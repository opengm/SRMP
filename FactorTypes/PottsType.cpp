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
#include "PottsType.h"

using namespace srmpLib;

PottsFactorType::PottsFactorType() 
	: buf(4096)
{
}

PottsFactorType::~PottsFactorType()
{
}

void PottsFactorType::InitFactor(Energy::NonSingletonFactor* A, double* user_data, unsigned flags)
{
	if (A->arity != 2 || A->nodes[0]->K != A->nodes[1]->K || *user_data < 0)
	{
		printf("Incorrect argument to PottsFactorType\n");
		exit(1);
	}
	A->data = (double*) buf.Alloc(sizeof(double));
	*A->data = *user_data;
}

double PottsFactorType::GetCost(Energy::NonSingletonFactor* A)
{
	return (A->nodes[1]->solution == A->nodes[0]->solution) ? 0 : (*A->data);
}

void PottsFactorType::InitEdge(Energy::Edge* e)
{
}

bool PottsFactorType::PrepareFactor(Energy::NonSingletonFactor* A)
{
	if (A->first_in) return false; // if there are incoming edges then convert factor to general 'PairwiseType'
	if (!A->first_out || !A->first_out->next_out || A->first_out->next_out->next_out
		|| A->first_out->B == A->first_out->next_out->B) return false; // degenerate cases, probably will never occur
	return true;
}

void PottsFactorType::ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* theta)
{
	// this function should never be called since incoming edges are not allowed (as specified in PrepareFactor())
	printf("Error: ComputePartialReparameterization() should not be called for this type. (Trying to save non-standard factor?");
	exit(1);
}


double PottsFactorType::SendMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	int k, K = e->B->K;
	double lambda = *A->data;

	Energy::Edge* e_rev = (e->next_out) ? e->next_out : A->first_out;
	assert (e_rev != e); // holds due to PrepareFactor()
	double* m_rev = e_rev->m;

	// need to set e->m[a] = min_{b} { cost(a,b) - m_rev[b] } and then normalize e->m
	// find the smallest value of (-m_rev[k])
	int k_best = 0;
	for (k=1; k<K; k++)
	{
		if (m_rev[k_best] < m_rev[k]) k_best = k;
	}
	double delta = -m_rev[k_best];
	e->m[k_best] = 0;
	for (k=0; k<K; k++)
	{
		if (k == k_best) continue;
		e->m[k] = -m_rev[k] - delta;
		if (e->m[k] > lambda) e->m[k] = lambda;
	}

	return delta;
}

void PottsFactorType::SendRestrictedMessage(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;
	int K = e->B->K, k_other = ((Energy::Node*)e->B == A->nodes[0]) ? A->nodes[1]->solution : A->nodes[0]->solution;
	double lambda = *A->data;

	memset(e->m, 0, K*sizeof(double));
	e->m[k_other] = -lambda;
}


double PottsFactorType::SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution)
{
	int K = A->nodes[0]->K;
	double* rep0 = A->nodes[0]->rep;
	double* rep1 = A->nodes[1]->rep;
	double lambda = *A->data;
	int weight0, weight1;
	if ((Energy::Node*)A->first_out->B == A->nodes[0])
		{ weight0 = A->first_out->weight_forward; weight1 = A->first_out->next_out->weight_forward; }
	else
		{ weight1 = A->first_out->weight_forward; weight0 = A->first_out->next_out->weight_forward; }
	

	int k, k0=0, k1=0, k01=0;
	double v0_min = rep0[0], v1_min = rep1[0], v01_min = rep0[0] + rep1[0];
	for (k=1; k<K; k++)
	{
		if (v0_min > rep0[k]) { v0_min = rep0[k]; k0 = k; }
		if (v1_min > rep1[k]) { v1_min = rep1[k]; k1 = k; }
		if (v01_min > rep0[k]+rep1[k]) { v01_min = rep0[k]+rep1[k]; k01 = k; }
	}
	double delta = v0_min + v1_min + lambda;
	if (delta > v01_min) delta = v01_min;

	if (set_solution)
	{
		if (A->nodes[0]->solution < 0)
		{
			if (A->nodes[1]->solution < 0)
			{
				if (delta < v01_min) { A->nodes[0]->solution = k0; A->nodes[1]->solution = k1; }
				else                 { A->nodes[0]->solution = A->nodes[1]->solution = k01; }
			}
			else
			{
				A->nodes[0]->solution = (rep0[A->nodes[1]->solution] < v0_min + lambda) ? A->nodes[1]->solution : k0;
			}
		}
		else if (A->nodes[1]->solution < 0)
		{
			A->nodes[1]->solution = (rep1[A->nodes[0]->solution] < v1_min + lambda) ? A->nodes[0]->solution : k1;
		}
	}

	int total_weight = A->weight_forward + weight0 + weight1;
	double total_weight_inv = 1.0 / total_weight;
	double rho0 = weight0*total_weight_inv, rho1 = weight1*total_weight_inv;

	for (k=0; k<K; k++)
	{
		double v, rep0_old = rep0[k], rep1_old = rep1[k];

		v = v1_min + lambda; if (v > rep1_old) v = rep1_old;
		rep0[k] = rho0*(rep0_old + v - delta);

		v = v0_min + lambda; if (v > rep0_old) v = rep0_old;
		rep1[k] = rho1*(rep1_old + v - delta);
	}

	return delta;
}

