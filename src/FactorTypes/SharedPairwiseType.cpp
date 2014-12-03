/*
    This file is part of SRMP software.

    SRMP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SRMP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SRMP.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <assert.h>
#include <srmp/FactorTypes/SharedPairwiseType.h>

using namespace srmpLib;

SharedPairwiseFactorType::SharedPairwiseFactorType(int _K1, int _K2, double* _costs) 
	: K1(_K1), K2(_K2),
	  buf(4096)
{
	costs = new double[K1*K2];
	memcpy(costs, _costs, K1*K2*sizeof(double));
}

SharedPairwiseFactorType::~SharedPairwiseFactorType()
{
	delete [] costs;
}

void SharedPairwiseFactorType::InitFactor(Energy::NonSingletonFactor* A, double* user_data)
{
	if (A->arity != 2 || A->nodes[0]->K != K1 || A->nodes[1]->K != K2 || *user_data <= 0)
	{
		printf("Incorrect argument to SharedPairwiseFactorType\n");
		exit(1);
	}
	A->data = (double*) buf.Alloc(sizeof(double));
	*A->data = *user_data;
}

double SharedPairwiseFactorType::GetCost(Energy::NonSingletonFactor* A)
{
	return (*A->data) * costs[A->nodes[1]->solution + A->nodes[1]->K*A->nodes[0]->solution];
}

void SharedPairwiseFactorType::InitEdge(Energy::Edge* e)
{
    (void)e; // fix unused variable compiler warning if build with NDEBUG defined

	//Energy::NonSingletonFactor* A = e->A;

	//assert (!A->first_out->next_out // single outgoing edge
	//	|| (!A->first_out->next_out->next_out && A->first_out->B != A->first_out->next_out->B) ); // two outgoing edges to different nodes
    assert (!e->A->first_out->next_out // single outgoing edge
        || (!e->A->first_out->next_out->next_out && e->A->first_out->B != e->A->first_out->next_out->B) ); // two outgoing edges to different nodes
}

bool SharedPairwiseFactorType::PrepareFactor(Energy::NonSingletonFactor*)
{
	return true;
}

void SharedPairwiseFactorType::ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* theta)
{
	Energy::Edge* e;
	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;

	double lambda = *A->data;
	for (x=0; x<A->K; x++) theta[x] = lambda*costs[x];

	for (e=A->first_out; e; e=e->next_out)
	{
		if ((Energy::Node*)e->B == A->nodes[0])
		{
			for (x=0; x<X; x++)
			for (y=0; y<Y; y++)
			{
				theta[x*Y + y] -= e->m[x];
			}
		}
		else
		{
			for (x=0; x<X; x++)
			for (y=0; y<Y; y++)
			{
				theta[x*Y + y] -= e->m[y];
			}
		}
	}
}


double SharedPairwiseFactorType::SendMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	Energy::Edge* e2;

	double lambda = *A->data;
	double lambda_inv = 1/lambda;

	double* theta;
	double* m_rev;

	if (!A->first_in) theta = costs;
	else
	{
		int a, KA = e->A->K;
		theta = (double*) rbuf.Alloc(KA*sizeof(double));
		memcpy(theta, costs, KA*sizeof(double));
		for (e2=A->first_in; e2; e2=e2->next_in)
		{
			for (a=0; a<KA; a++) theta[a] += lambda_inv * e2->m[a];
		}
	}

	int K_rev = K1 + K2 - e->B->K;
	m_rev = (double*) rbuf2.Alloc(K_rev*sizeof(double));
	e2 = (e->next_out) ? e->next_out : A->first_out;
	if (e2 != e)
	{
		int x;
		for (x=0; x<K_rev; x++) m_rev[x] = lambda_inv * e2->m[x];
	}
	else memset(m_rev, 0, K_rev*sizeof(double));

	int x, y, X = K1, Y = K2;
	if ((Energy::Node*)e->B == A->nodes[0])
	{
		for (x=0; x<X; x++)
		{
/*
			double v = theta[x*Y + 0] - m_rev[0];
			for (y=1; y<Y; y++)
			{
				if (v > theta[x*Y + y] - m_rev[y]) v = theta[x*Y + y] - m_rev[y];
			}
*/
			double* ptr = theta + x*Y;
			double* m_rev_ptr = m_rev;
			double v = (*ptr) - (*m_rev_ptr);
			for (ptr++, m_rev_ptr++; m_rev_ptr<m_rev+Y; ptr++, m_rev_ptr++)
			{
				if (v > (*ptr) - (*m_rev_ptr)) v = (*ptr) - (*m_rev_ptr);
			}

			e->m[x] = lambda * v;
		}
	}
	else
	{
		for (y=0; y<Y; y++)
		{
/*
			double v = theta[0*Y + y] - m_rev[0];
			for (x=1; x<X; x++)
			{
				if (v > theta[x*Y + y] - m_rev[x]) v = theta[x*Y + y] - m_rev[x];
			}
*/
			double* ptr = theta + y;
			double* m_rev_ptr = m_rev;
			double v = (*ptr) - (*m_rev_ptr);
			for (ptr+=Y, m_rev_ptr++; m_rev_ptr<m_rev+X; ptr+=Y, m_rev_ptr++)
			{
				if (v > (*ptr) - (*m_rev_ptr)) v = (*ptr) - (*m_rev_ptr);
			}

			e->m[y] = lambda * v;
		}
	}
	double delta = e->m[0];
	for (x=1; x<e->B->K; x++)
	{
		if (delta > e->m[x]) delta = e->m[x];
	}
	for (x=0; x<e->B->K; x++) e->m[x] -= delta;

	return delta;
}

void SharedPairwiseFactorType::SendRestrictedMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	Energy::Edge* e2;
	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;
	double lambda = *A->data;

	if ((Energy::Node*)e->B == A->nodes[0])
	{
		assert(A->nodes[0]->solution < 0 && A->nodes[1]->solution >= 0);
		y = A->nodes[1]->solution;
		for (x=0; x<X; x++)
		{
			e->m[x] = lambda * costs[x*Y + y];
		}
		for (e2=A->first_in; e2; e2=e2->next_in)
		{
			for (x=0; x<X; x++)
			{
				e->m[x] += e2->m[x*Y + y];
			}
		}
	}
	else
	{
		assert(A->nodes[0]->solution >= 0 && A->nodes[1]->solution < 0);
		x = A->nodes[0]->solution;
		for (y=0; y<Y; y++)
		{
			e->m[y] = lambda * costs[x*Y + y];
		}
		for (e2=A->first_in; e2; e2=e2->next_in)
		{
			for (y=0; y<Y; y++)
			{
				e->m[y] += e2->m[x*Y + y];
			}
		}
	}
}



double SharedPairwiseFactorType::SendMPLPMessages(Energy::NonSingletonFactor*, bool)
{
	printf("SharedPairwiseFactorType::SendMPLPMessages() is not implemented\n");
	exit(1);
	return 0;
}
