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
#include <string.h>
#include <assert.h>
#include "PatternType.h"

using namespace srmpLib;

inline int CAST(void* p)
{
	union { int i; void* p;	} t;
	t.i = 0;
	t.p = p;
	return t.i;
}

inline void* CAST(int i)
{
	union { int i; void* p;	} t;
	t.p = NULL;
	t.i = i;
	return t.p;
}

//////////////////////////////////////////////////////////////////////////////////

PatternFactorType::PatternFactorType() 
	: buf(4096)
{
}

PatternFactorType::~PatternFactorType()
{
}

void PatternFactorType::InitFactor(Energy::NonSingletonFactor* A, double* user_data, unsigned flags)
{
	Input* input = (Input*) user_data;
	if (input->cost > 0)
	{
		printf("Incorrect argument to PatternFactorType\n");
		exit(1);
	}
	FactorData* D = (FactorData*) buf.Alloc(sizeof(FactorData));
	A->data = (double*) D;
	D->cost = input->cost;
	if (flags & FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY)
	{
		D->pattern = input->pattern;
	}
	else
	{
		D->pattern = (int*) buf.Alloc(A->arity*sizeof(int));
		memcpy(D->pattern, input->pattern, A->arity*sizeof(int));
	}
}

double PatternFactorType::GetCost(Energy::NonSingletonFactor* A)
{
	FactorData* D = (FactorData*) A->data;
	int i;
	for (i=0; i<A->arity; i++)
	{
		if (A->nodes[i]->solution != D->pattern[i]) return 0;
	}
	return D->cost;
}

void PatternFactorType::InitEdge(Energy::Edge* e)
{
}

bool PatternFactorType::PrepareFactor(Energy::NonSingletonFactor* A)
{
	if (A->first_in) return false; // if there are incoming edges then convert factor to general 'PairwiseType'

	int i, n = A->arity;
	Energy::Edge* e;

	for (e=A->first_out; e; e=e->next_out)
	{
		if (e->B->arity != 1) return false;
	}

	// set e->send_message_data so that D->pattern[(int)e->send_message_data] == e->B
	// temporarily use Node::arity for this purpose
	bool degenerate_case = false;
	for (i=0; i<n; i++) A->nodes[i]->arity = i;
	for (e=A->first_out; e; e=e->next_out)
	{
		if (e->B->arity < 0) { degenerate_case = true; break; }
		e->send_message_data = CAST(e->B->arity);
		e->B->arity = -1;
	}
	for (i=0; i<n; i++)
	{
		if (A->nodes[i]->arity >= 0) { degenerate_case = true; break; }
	}
	for (i=0; i<n; i++) A->nodes[i]->arity = 1;

	FactorData* D = (FactorData*) A->data;
	D->counter = 0;
	D->e_last = A->first_out;

	return true;
}

void PatternFactorType::RecomputeFactorData(Energy::NonSingletonFactor* A)
{
	FactorData* D = (FactorData*) A->data;
	int i, k, n = A->arity;
	Energy::Edge* e;

	D->sum_min = 0;
	D->sum_pattern = D->cost;
	for (e=A->first_out; e; e=e->next_out)
	{
		if (e == D->e_last) continue;
		i = CAST(e->send_message_data);

		double v_min = -e->m[0];
		for (k=1; k<e->B->K; k++)
		{
			if (v_min > -e->m[k]) v_min = -e->m[k];
		}
		D->sum_min += v_min;
		D->sum_pattern += -e->m[D->pattern[i]];
	}
}

void PatternFactorType::ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* theta)
{
	// this function should never be called since incoming edges are not allowed (as specified in PrepareFactor())
	printf("Error: ComputePartialReparameterization() should not be called for this type. (Trying to save non-standard factor?");
	exit(1);
}


double PatternFactorType::SendMessage(Energy::Edge* e0)
{
	Energy::NonSingletonFactor* A = e0->A;
	FactorData* D = (FactorData*) A->data;
	int i, k, n = A->arity;

	if (D->counter == 0)
	{
		RecomputeFactorData(A);
		D->counter = 10*n;
	}
	else D->counter --;

	// change D->e_last to e0, update D->sum_min and D->sum_pattern accordingly
	// first, add values for D->e_last
	double v_min = -D->e_last->m[0];
	for (k=1; k<D->e_last->B->K; k++)
	{
		if (v_min > -D->e_last->m[k]) v_min = -D->e_last->m[k];
	}
	D->sum_min += v_min;
	i = CAST(D->e_last->send_message_data);
	D->sum_pattern += -D->e_last->m[D->pattern[i]];
	D->e_last = e0;

	// now subtract values for e0
	v_min = -e0->m[0];
	for (k=1; k<e0->B->K; k++)
	{
		if (v_min > -e0->m[k]) v_min = -e0->m[k];
	}
	D->sum_min -= v_min;
	i = CAST(e0->send_message_data);
	D->sum_pattern -= -e0->m[D->pattern[i]];

	// now main computations
	if (D->sum_min < D->sum_pattern)
	{
		for (k=0; k<e0->B->K; k++) e0->m[k] = 0;
		return D->sum_min;
	}
	else
	{
		for (k=0; k<e0->B->K; k++) e0->m[k] = D->sum_min - D->sum_pattern;
		e0->m[D->pattern[i]] = 0;
		return D->sum_pattern;
	}
}

void PatternFactorType::SendRestrictedMessage(Energy::Edge* e0)
{
	Energy::NonSingletonFactor* A = e0->A;
	FactorData* D = (FactorData*) A->data;
	int i, k, n = A->arity;
	Energy::Edge* e;

	for (k=0; k<e0->B->K; k++) e0->m[k] = 0;

	// check whether labeled nodes are consistent with 'pattern'
	for (i=0; i<n; i++)
	{
		if (A->nodes[i]->solution >= 0 && A->nodes[i]->solution != D->pattern[i])
		{
			return; // not consistent!
		}
	}

	double sum_min = 0, sum_pattern = D->cost; // sums over unlabeled nodes (and also excluding e0->B)
	for (e=A->first_out; e; e=e->next_out)
	{
		if (((Energy::Node*)e->B)->solution >= 0 || e == e0) continue;
		i = CAST(e->send_message_data);

		double v_min = -e->m[0];
		for (k=1; k<e->B->K; k++)
		{
			if (v_min > -e->m[k]) v_min = -e->m[k];
		}
		sum_min += v_min;
		sum_pattern += -e->m[D->pattern[i]];
	}

	if (sum_pattern < sum_min)
	{
		i = CAST(e0->send_message_data);
		e0->m[D->pattern[i]] = sum_pattern - sum_min;
	}
}


double PatternFactorType::SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution)
{
	FactorData* D = (FactorData*) A->data;
	int _i, k, n = A->arity;
	Energy::Node* i;
	Energy::Edge* e;

	double* v_min_array = (double*) rbuf.Alloc(n*sizeof(double));
	double sum_min = 0, sum_pattern = D->cost;

	for (_i=0; _i<n; _i++)
	{
		i = A->nodes[_i];
		double v_min = i->rep[0];
		for (k=1; k<i->K; k++)
		{
			if (v_min > i->rep[k]) v_min = i->rep[k];
		}
		v_min_array[_i] = v_min;
		sum_min += v_min;
		sum_pattern += i->rep[D->pattern[_i]];
	}

	if (set_solution)
	{
		_i = 0;
		if (sum_min >= sum_pattern)
		{
			for (_i=0; _i<n; _i++)
			{
				i = A->nodes[_i];
				if (i->solution >= 0 && i->solution != D->pattern[_i]) break;
			}
		}
		if (_i==n) // pattern cost is smaller
		{
			for (_i=0; _i<n; _i++)
			{
				i = A->nodes[_i];
				if (i->solution < 0) i->solution = D->pattern[_i];
			}
		}
		else // default cost is smaller
		{
			for (_i=0; _i<n; _i++)
			{
				i = A->nodes[_i];
				if (i->solution >= 0) continue;
				for (k=0; ; k++)
				{
					if (i->rep[k] == v_min_array[_i]) break;
				}
				i->solution = k;
			}
		}
	}

	double delta = (sum_min < sum_pattern) ? sum_min : sum_pattern;
	sum_min -= delta;
	sum_pattern -= delta;

	int w = A->weight_forward;
	for (e=A->first_out; e; e=e->next_out) w += e->weight_forward;
	double w_inv = 1.0 / w;
	for (e=A->first_out; e; e=e->next_out)
	{
		double rho = e->weight_forward * w_inv;
		_i = CAST(e->send_message_data);
		i = A->nodes[_i];
		double sum_min2 = sum_min - v_min_array[_i];

		for (k=0; k<i->K; k++)
		{
			if (k == D->pattern[_i])
			{
				if (sum_pattern < i->rep[k] + sum_min2) i->rep[k] = rho*sum_pattern;
				else                                    i->rep[k] = rho*(i->rep[k] + sum_min2);
			}
			else
			{
				i->rep[k] = rho*(i->rep[k] + sum_min2);
			}
		}
	}

	return delta;
}

