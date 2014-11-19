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
#include "PairwiseType.h"

#define SRMP_WHICH_NODE user1

PairwiseFactorType::PairwiseFactorType() : buf(4096)
{
}

PairwiseFactorType::~PairwiseFactorType()
{
}

void PairwiseFactorType::InitFactor(Energy::NonSingletonFactor* A, double* user_data, unsigned flags)
{
	if ((flags & FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY) && user_data)
	{
		A->data = user_data;
	}
	else
	{
		A->data = (double*) buf.Alloc(A->K*sizeof(double));
		if (user_data) memcpy(A->data, user_data, A->K*sizeof(double));
		else           memset(A->data, 0, A->K*sizeof(double));
	}
}

double PairwiseFactorType::GetCost(Energy::NonSingletonFactor* A)
{
	if (!A->data) return 0;
	return A->data[A->nodes[1]->solution + A->nodes[1]->K*A->nodes[0]->solution];
}

void PairwiseFactorType::InitEdge(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;

	assert (!A->first_out->next_out // single outgoing edge
		|| (!A->first_out->next_out->next_out && A->first_out->B != A->first_out->next_out->B) ); // two outgoing edges to different nodes

	if ((Energy::Node*)e->B == A->nodes[0]) e->SRMP_WHICH_NODE = 0;
	else                                    e->SRMP_WHICH_NODE = 1;
}

bool PairwiseFactorType::PrepareFactor(Energy::NonSingletonFactor* A)
{
	return true;
}

void PairwiseFactorType::ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* theta)
{
	Energy::Edge* e;
	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;

	memcpy(theta, A->data, A->K*sizeof(double));
	for (e=A->first_out; e; e=e->next_out)
	{
		if (e->SRMP_WHICH_NODE == 0)
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


double PairwiseFactorType::SendMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	Energy::Edge* e2;

	double* theta;
	double* m_rev;

	if (!A->first_in) theta = A->data;
	else
	{
		int a, KA = A->K;
		theta = (double*) rbuf.Alloc(KA*sizeof(double));
		memcpy(theta, A->data, KA*sizeof(double));
		for (e2=A->first_in; e2; e2=e2->next_in)
		{
			for (a=0; a<KA; a++) theta[a] += e2->m[a];
		}
	}

	e2 = (e->next_out) ? e->next_out : A->first_out;
	if (e2 != e) m_rev = e2->m;
	else
	{
		int K_rev = A->nodes[1 - e->SRMP_WHICH_NODE]->K;
		m_rev = (double*) rbuf2.Alloc(K_rev*sizeof(double));
		memset(m_rev, 0, K_rev*sizeof(double));
	}

	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;
	if (e->SRMP_WHICH_NODE == 0)
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

			e->m[x] = v;
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

			e->m[y] = v;
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

void PairwiseFactorType::SendRestrictedMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	Energy::Edge* e2;
	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;

	assert(A->nodes[e->SRMP_WHICH_NODE]->solution < 0 && A->nodes[1 - e->SRMP_WHICH_NODE]->solution >= 0);

	if (e->SRMP_WHICH_NODE == 0)
	{
		y = A->nodes[1]->solution;
		for (x=0; x<X; x++)
		{
			e->m[x] = A->data[x*Y + y];
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
		x = A->nodes[0]->solution;
		for (y=0; y<Y; y++)
		{
			e->m[y] = A->data[x*Y + y];
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


double PairwiseFactorType::_SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution)
{
	int total_weight = A->weight_forward;
	int a, KA = A->K;
	double delta = 0;

	Energy::Edge* e;

	double* theta = (double*) rbuf.Alloc(KA*sizeof(double) + 4*A->arity*sizeof(int));
	memcpy(theta, A->data, KA*sizeof(double));
	for (e=A->first_in; e; e=e->next_in)
	{
		for (a=0; a<KA; a++) theta[a] += e->m[a];
	}

	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;
	for (e=A->first_out; e; e=e->next_out)
	{
		if (e->SRMP_WHICH_NODE == 0)
		{
			for (x=0, a=0; x<X; x++)
			for (y=0; y<Y; y++, a++)
			{
				theta[a] += e->B->rep[x];
			}
		}
		else
		{
			for (x=0, a=0; x<X; x++)
			for (y=0; y<Y; y++, a++)
			{
				theta[a] += e->B->rep[y];
			}
		}
		total_weight += e->weight_forward;
	}

	if (set_solution) A->ComputeRestrictedMinimum(theta, (int*)(theta+KA));

	delta = theta[0];
	for (a=1; a<KA; a++)
	{
		if (delta > theta[a]) delta = theta[a];
	}

	double total_weight_inv = 1.0 / total_weight;

	if (A->rep) memcpy(A->rep, theta, A->K*sizeof(double));

	for (e=A->first_out; e; e=e->next_out)
	{
		double rho = e->weight_forward * total_weight_inv;

		if (e->SRMP_WHICH_NODE == 0)
		{
			for (x=0; x<X; x++)
			{
/*				double v = theta[x*Y + 0];
				for (y=1; y<Y; y++)
				{
					if (v > theta[x*Y + y]) v = theta[x*Y + y];
				}
*/
				double* ptr = theta + x*Y;
				double* ptr_end = ptr + Y;
				double v = (*ptr);
				for (ptr++; ptr<ptr_end; ptr++)
				{
					if (v > (*ptr)) v = (*ptr);
				}

				e->B->rep[x] = rho * (v - delta);
			}
			if (A->rep)
			{
				for (a=0, x=0; x<X; x++)
				for (y=0; y<Y; y++, a++)
				{
					A->rep[a] -= e->B->rep[x];
				}
			}
		}
		else
		{
			for (y=0; y<Y; y++)
			{
/*
				double v = theta[0*Y + y];
				for (x=1; x<X; x++)
				{
					if (v > theta[x*Y + y]) v = theta[x*Y + y];
				}
*/
				double* ptr = theta + y;
				double v = (*ptr);
				for (ptr+=Y; ptr<theta+KA; ptr+=Y)
				{
					if (v > (*ptr)) v = (*ptr);
				}

				e->B->rep[y] = rho * (v - delta);
			}
			if (A->rep)
			{
				for (a=0, x=0; x<X; x++)
				for (y=0; y<Y; y++, a++)
				{
					A->rep[a] -= e->B->rep[y];
				}
			}
		}
	}

	return delta;
}

double PairwiseFactorType::SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution)
{
	if (set_solution || A->rep || !A->first_out || !A->first_out->next_out || A->first_out->weight_forward==0 || A->first_out->next_out->weight_forward==0)
	{
		return PairwiseFactorType::_SendMPLPMessages(A, set_solution);
	}
	assert (!A->first_in); // since A->rep is NULL

	Energy::Edge* e;

	double* theta = A->data;
	double delta = 0;

	double* rep[2];
	int weight[2];
	for (e=A->first_out; e; e=e->next_out)
	{
		rep[e->SRMP_WHICH_NODE] = e->B->rep;
		weight[e->SRMP_WHICH_NODE] = e->weight_forward;
	}

	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;

	double* rep0_old = (double*) rbuf.Alloc(X*sizeof(double));
	memcpy(rep0_old, rep[0], X*sizeof(double));

	int total_weight = A->weight_forward + weight[0] + weight[1];
	double total_weight_inv = 1.0 / total_weight;
	double rho[2] = { weight[0] * total_weight_inv, weight[1] * total_weight_inv };

	for (x=0; x<X; x++)
	{
/*
		double v = theta[x*Y + 0] + rep[1][0];
		for (y=1; y<Y; y++)
		{
			if (v > theta[x*Y + y] + rep[1][y]) v = theta[x*Y + y] + rep[1][y];
		}
*/
		double* ptr = theta + x*Y;
		double* ptr1 = rep[1];
		double v = (*ptr) + (*ptr1);
		for (ptr++, ptr1++; ptr1<rep[1]+Y; ptr++, ptr1++)
		{
			if (v > (*ptr) + (*ptr1)) v = (*ptr) + (*ptr1);
		}

		v += rep[0][x];
		if (x == 0 || delta > v) delta = v;
		rep[0][x] = rho[0] * v;
	}
	for (y=0; y<Y; y++)
	{
/*
		double v = theta[0*Y + y] + rep0_old[0];
		for (x=1; x<X; x++)
		{
			if (v > theta[x*Y + y] + rep0_old[x]) v = theta[x*Y + y] + rep0_old[x];
		}
*/
		double* ptr = theta + y;
		double* ptr0 = rep0_old;
		double v = (*ptr) + (*ptr0);
		for (ptr+=Y, ptr0++; ptr0<rep0_old+X; ptr+=Y, ptr0++)
		{
			if (v > (*ptr) + (*ptr0)) v = (*ptr) + (*ptr0);
		}

		v += rep[1][y];
		rep[1][y] = rho[1] * (v - delta);
	}

	for (x=0; x<X; x++) rep[0][x] -= rho[0] * delta;

	return delta;
}








#ifdef FHLAFA
double PairwiseFactorType::SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution)
{
	if (A->rep || !A->first_out || !A->first_out->next_out || A->first_out->weight_forward==0 || A->first_out->weight_forward==0)
	{
	}

	Energy::Edge* e;

	double* theta;
	double delta = 0;

	if (!A->first_in) theta = A->data;
	else
	{
		int a, KA = A->K;
		theta = (double*) rbuf.Alloc(KA*sizeof(double));
		memcpy(theta, A->data, KA*sizeof(double));
		for (e=A->first_in; e; e=e->next_in)
		{
			for (a=0; a<KA; a++) theta[a] += e->m[a];
		}
	}

	e = A->first_out;
	if (!e)
	{
		int a, KA = A->K;
		delta = theta[0];
		for (a=1; a<KA; a++)
		{
			if (delta > theta[a]) delta = theta[a];
		}
		if (A->rep) memcpy(A->rep, theta, KA*sizeof(double));
		return delta;
	}

	double* rep[2];
	int weight[2];
	rep[e->SRMP_WHICH_NODE] = e->B->rep;
	weight[e->SRMP_WHICH_NODE] = e->weight_forward;
	if (e->next_out)
	{
		e = e->next_out;
		rep[e->SRMP_WHICH_NODE] = e->B->rep;
		weight[e->SRMP_WHICH_NODE] = e->weight_forward;
	}
	else
	{
		int K = A->nodes[1 - e->SRMP_WHICH_NODE]->K;
		rep[1 - e->SRMP_WHICH_NODE] = (double*) rbuf2.Alloc(K*sizeof(double));
		memset(rep[1 - e->SRMP_WHICH_NODE], 0, K*sizeof(double));
		weight[1 - e->SRMP_WHICH_NODE] = 0;
	}

	int x, y, X = A->nodes[0]->K, Y = A->nodes[1]->K;

	if (A->rep)
	{
		int a, KA = A->K;
		memcpy(A->rep, theta, KA*sizeof(double));
		for (x=0, a=0; x<X; x++)
		for (y=0; y<Y; y++, a++)
		{
			A->rep[a] += rep[0][x] + rep[1][y];
		}
	}

	double* rep0_old = (double*) rbuf3.Alloc(X*sizeof(double));
	memcpy(rep0_old, rep[0], X*sizeof(double));

	int total_weight = A->weight_forward + weight[0] + weight[1];
	double total_weight_inv = 1.0 / total_weight;
	double rho[2] = { weight[0] * total_weight_inv, weight[1] * total_weight_inv };

	for (x=0; x<X; x++)
	{
/*
		double v = theta[x*Y + 0] + rep[1][0];
		for (y=1; y<Y; y++)
		{
			if (v > theta[x*Y + y] + rep[1][y]) v = theta[x*Y + y] + rep[1][y];
		}
*/
		double* ptr = theta + x*Y;
		double* ptr1 = rep[1];
		double v = (*ptr) + (*ptr1);
		for (ptr++, ptr1++; ptr1<rep[1]+Y; ptr++, ptr1++)
		{
			if (v > (*ptr) + (*ptr1)) v = (*ptr) + (*ptr1);
		}

		v += rep[0][x];
		if (x == 0 || delta > v) delta = v;
		rep[0][x] = rho[0] * v;
	}
	for (y=0; y<Y; y++)
	{
/*
		double v = theta[0*Y + y] + rep0_old[0];
		for (x=1; x<X; x++)
		{
			if (v > theta[x*Y + y] + rep0_old[x]) v = theta[x*Y + y] + rep0_old[x];
		}
*/
		double* ptr = theta + y;
		double* ptr0 = rep0_old;
		double v = (*ptr) + (*ptr0);
		for (ptr+=Y, ptr0++; ptr0<rep0_old+X; ptr+=Y, ptr0++)
		{
			if (v > (*ptr) + (*ptr0)) v = (*ptr) + (*ptr0);
		}

		v += rep[1][y];
		rep[1][y] = rho[1] * (v - delta);
	}

	if (A->rep)
	{
		int a, KA = A->K;
		for (x=0, a=0; x<X; x++)
		for (y=0; y<Y; y++, a++)
		{
			A->rep[a] -= rep[0][x] + rep[1][y] - delta*rho[0];
		}
	}

	for (x=0; x<X; x++) rep[0][x] -= rho[0] * delta;

	return delta;
}
#endif
