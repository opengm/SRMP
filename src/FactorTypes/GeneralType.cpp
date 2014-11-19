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
#include <srmp/FactorTypes/GeneralType.h>

using namespace srmpLib;

GeneralFactorType::GeneralFactorType() 
	: buf(4096)
{
}

GeneralFactorType::~GeneralFactorType()
{
}

void GeneralFactorType::InitFactor(Energy::NonSingletonFactor* A, double* user_data, unsigned flags)
{
	if (flags & FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY)
	{
		A->data = user_data;
	}
	else if (user_data)
	{
		A->data = (double*) buf.Alloc(A->K*sizeof(double));
		memcpy(A->data, user_data, A->K*sizeof(double));
	}
	else
	{
		A->data = NULL;
	}
}

double GeneralFactorType::GetCost(Energy::NonSingletonFactor* A)
{
	if (!A->data) return 0;

	int i, a = 0, K_factor=1;

	for (i=A->arity-1; ; i--)
	{
		a += K_factor*A->nodes[i]->solution;
		if (i==0) return A->data[a];
		K_factor *= A->nodes[i]->K;
	}
}

void SetEdgeTable(int nA, Energy::Node** A, int nB, Energy::Node** B, int* table, int* _buf)
{
	int i, j, b=0, k=0;
	int* K_array = _buf;
	int* labeling = _buf + nB;
	for (i=0; i<nB; i++)
	{
		K_array[i] = 1;
		for (j=nA-1; A[j] != B[i]; j--) K_array[i] *= A[j]->K;
	}
	memset(labeling, 0, nB*sizeof(int));
	table[0] = 0;
	while ( 1 )
	{
		for (i=nB-1; i>=0; i--)
		{
			if (labeling[i] < B[i]->K-1) break;
			k -= labeling[i]*K_array[i];
			labeling[i] = 0;
		}
		if (i<0) break;
		labeling[i] ++;
		k += K_array[i];
		table[++b] = k;
	}
}

void GeneralFactorType::InitEdge(Energy::Edge* e)
{
	int i, j;
	Energy::NonSingletonFactor* A = e->A;
	Energy::Factor* B = e->B;
	int* _buf = (int*) rbuf.Alloc(2*A->arity*sizeof(int) + A->arity*sizeof(Energy::Node*));

	Energy::Node** Anodes_sorted = A->nodes + A->is_unsorted*A->arity;
	Energy::Node** Cnodes = (Energy::Node**) (_buf + 2*A->arity);
	Energy::Node** Bnodes;
	Energy::Node** Bnodes_sorted;
	if (B->arity == 1) Bnodes = Bnodes_sorted = (Energy::Node**)&B;
	else { Bnodes = ((Energy::NonSingletonFactor*)B)->nodes; Bnodes_sorted = Bnodes + B->arity*B->is_unsorted; }

	e->send_message_data = buf.Alloc((B->K + A->K/B->K)*sizeof(int));
	int* TB = (int*) e->send_message_data;
	int* TC = TB + B->K;

	SetEdgeTable(A->arity, A->nodes, B->arity, Bnodes, TB, _buf);

	for (i=j=0; ; i++)
	{
		while (j < B->arity && Anodes_sorted[i] == Bnodes_sorted[j]) { i++; j++; }
		if (j == B->arity) break;
		Cnodes[i-j] = Anodes_sorted[i];
	}
	for ( ; i<A->arity; i++) Cnodes[i-j] = Anodes_sorted[i];

	SetEdgeTable(A->arity, A->nodes, i-j, Cnodes, TC, _buf);
}

bool GeneralFactorType::PrepareFactor(Energy::NonSingletonFactor* A)
{
	return true;
}

void GeneralFactorType::ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* theta)
{
	Energy::Edge* e;
	int b, c, KA = A->K;

	if (A->data) memcpy(theta, A->data, KA*sizeof(double));
	else         memset(theta, 0, KA*sizeof(double));
	for (e=A->first_out; e; e=e->next_out)
	{
		int KB = e->B->K;
		int KC = KA / KB;
		int* TB = (int*) e->send_message_data;
		int* TC = TB + KB;
		for (b=0; b<KB; b++)
		for (c=0; c<KC; c++)
		{
			theta[TB[b] + TC[c]] -= e->m[b];
		}
	}
}

double GeneralFactorType::SendMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	int a, b, c, KA = A->K;
	double delta = 0;

	Energy::Edge* e2;

	double* theta = (double*) rbuf.Alloc(KA*sizeof(double));
	if (A->data) memcpy(theta, A->data, KA*sizeof(double));
	else         memset(theta, 0, KA*sizeof(double));
	for (e2=A->first_in; e2; e2=e2->next_in)
	{
		for (a=0; a<KA; a++) theta[a] += e2->m[a];
	}
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		if (e2 == e) continue;
		int KB = e2->B->K;
		int KC = KA / KB;
		int* TB = (int*) e2->send_message_data;
		int* TC = TB + KB;
		for (b=0; b<KB; b++)
		for (c=0; c<KC; c++)
		{
			theta[TB[b] + TC[c]] -= e2->m[b];
		}
	}
	int KB = e->B->K;
	int KC = KA / KB;
	int* TB = (int*) e->send_message_data;
	int* TC = TB + KB;
	for (b=0; b<KB; b++)
	{
		double v_min = theta[TB[b]]; // TC[c] == 0
		for (c=1; c<KC; c++)
		{
			if (v_min > theta[TB[b] + TC[c]]) v_min = theta[TB[b] + TC[c]];
		}
		e->m[b] = v_min;
		if (b==0 || delta>v_min) delta = v_min;
	}
	for (b=0; b<KB; b++) e->m[b] -= delta;

	return delta;
}

void GeneralFactorType::SendRestrictedMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	int a, b, c, KA = A->K;

	Energy::Edge* e2;

	double* theta = (double*) rbuf.Alloc(KA*sizeof(double) + 4*A->arity*sizeof(int));
	int* _buf = (int*)(theta + KA);
	if (A->data) memcpy(theta, A->data, KA*sizeof(double));
	else         memset(theta, 0, KA*sizeof(double));
	for (e2=A->first_in; e2; e2=e2->next_in)
	{
		for (a=0; a<KA; a++) theta[a] += e2->m[a];
	}
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		if (e2 == e) continue;
		int KB = e2->B->K;
		int KC = KA / KB;
		int* TB = (int*) e2->send_message_data;
		int* TC = TB + KB;
		for (b=0; b<KB; b++)
		for (c=0; c<KC; c++)
		{
			theta[TB[b] + TC[c]] -= e2->m[b];
		}
	}

	/////////////////////////////////////////////////////////////

	A->ComputeRestrictedMinimum(theta, e->B, e->m, _buf);
}





double GeneralFactorType::SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution)
{
	int total_weight = A->weight_forward;
	int a, b, c, KA = A->K;
	double delta = 0;

	Energy::Edge* e;

	double* theta = (double*) rbuf.Alloc(KA*sizeof(double) + 4*A->arity*sizeof(int));
	if (A->data) memcpy(theta, A->data, KA*sizeof(double));
	else         memset(theta, 0, KA*sizeof(double));
	for (e=A->first_in; e; e=e->next_in)
	{
		for (a=0; a<KA; a++) theta[a] += e->m[a];
	}
	for (e=A->first_out; e; e=e->next_out)
	{
		int KB = e->B->K;
		int KC = KA / KB;
		int* TB = (int*) e->send_message_data;
		int* TC = TB + KB;
		for (b=0; b<KB; b++)
		for (c=0; c<KC; c++)
		{
			theta[TB[b] + TC[c]] += e->B->rep[b];
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
		int KB = e->B->K;
		int KC = KA / KB;
		int* TB = (int*) e->send_message_data;
		int* TC = TB + KB;
		for (b=0; b<KB; b++)
		{
			double v_min = theta[TB[b]]; // TC[c] == 0
			for (c=1; c<KC; c++)
			{
				if (v_min > theta[TB[b] + TC[c]]) v_min = theta[TB[b] + TC[c]];
			}
			e->B->rep[b] = rho * (v_min - delta);
		}
		if (A->rep)
		{
			for (b=0; b<KB; b++)
			{
				double v = e->B->rep[b];
				for (c=0; c<KC; c++)
				{
					A->rep[TB[b] + TC[c]] -= v;
				}
			}
		}
	}

	return delta;
}
	
