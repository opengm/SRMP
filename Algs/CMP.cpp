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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "util.h"
#include "../SRMP.h"




double Energy::InitCMP(Sequence& seq, Options& options)
{
	int t;
	Node* i;
	Factor* A;
	NonSingletonFactor* _A;
	Edge* e;

	// set 'separator_num' and 'order'
	int separator_num = 0;
	double LB_init = 0;
	for (i=nodes; i<nodes+node_num; i++) { if (!i->is_removed) separator_num ++; }
	for (_A=factors.ScanFirst(); _A; _A=factors.ScanNext())
	{
		if (_A->is_removed) continue;
		_A->tmp1 = 1;
		if (_A->first_in) separator_num ++;
		else if (!_A->first_out)
		{
			_A->rep = NULL;
			LB_init += SRMP_SEND_MPLP_MESSAGES(_A, false);
		}
	}

	seq.num = separator_num;
	seq.arr = new FactorPtr[separator_num];

	separator_num = 0;
	for (i=nodes; i<nodes+node_num; i++) { if (!i->is_removed) seq.arr[separator_num++].A = i; }
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (!A->is_removed && A->first_in) seq.arr[separator_num++].A = A;
	}

	SortSequence(seq, options);

	seq.arity_max = 0;
	seq.K_max = 0;
	for (t=separator_num-1; t>=0; t--)
	{
		A = seq.arr[t].A;

		if (seq.arity_max < A->arity) seq.arity_max = A->arity;
		if (seq.K_max < A->K) seq.K_max = A->K;

		A->compute_bound = (A->arity==1) ? 1 : A->tmp1; // if A->arity==1 then A->tmp1 hasn't been set
		A->tmp1 = 0;
		for (e=A->first_in; e; e=e->next_in)
		{
			e->compute_bound = e->A->tmp1;
			e->A->tmp1 = 0;
		}
	}

	// set weights
	for (t=0; t<separator_num; t++)
	{
		A = seq.arr[t].A;
		for (e=A->first_in; e; e=e->next_in)
		{
			e->weight_forward = 1;
		}
		A->weight_forward = 1;
	}

	return LB_init;
}

double Energy::SolveCMP(Options& options)
{
	double time_start = get_time();

	Sequence seq;
	double LB_init = InitCMP(seq, options);

	Node* i;
	int t, iter, iter_max = options.iter_max;
	Factor* B;
	Edge* e;

	double v, LB = 0, LB_prev = 0;
	double cost_fw;

	double* theta = new double[seq.K_max];
	char* restricted_buf = new char[2*seq.K_max*sizeof(double) + 4*seq.arity_max*sizeof(int)];

	for (iter=0; iter<iter_max; iter++)
	{
		if ((iter > 1 && LB<LB_prev + options.eps) || get_time()-time_start > options.time_max)
			iter_max = iter+1; // last iteration, with compute_solution=true
		LB_prev = LB;

		bool compute_solution = false;
		if (iter != 0 && iter % options.compute_solution_period == 0) compute_solution = true;
		if (iter == iter_max-1) compute_solution = true;


		// forward pass
		if (compute_solution)
		{
			for (i=nodes; i<nodes+node_num; i++) i->solution = -1;
		}
		LB = LB_init;
		for (t=0; t<seq.num; t++)
		{
			B = seq.arr[t].A;
			int b, K = B->K;

			if (B->arity == 1) memcpy(theta, B->data, K*sizeof(double));
			else SRMP_COMPUTE_PARTIAL_REPARAMETERIZATION((NonSingletonFactor*)B, theta);

			int w_total = B->weight_forward;
			int w_compute_bound = (B->compute_bound) ? B->weight_forward : 0;
			for (e=B->first_in; e; e=e->next_in)
			{
				v = SRMP_SEND_MESSAGE(e);
				if (e->compute_bound) { LB += v; w_compute_bound += e->weight_forward; }
				for (b=0; b<K; b++) theta[b] += e->m[b]; 
				w_total += e->weight_forward;
			}

			if (compute_solution)
			{
				ComputeSolution(B, restricted_buf);
			}

			double p = 1.0 / w_total;
			for (b=0; b<K; b++) theta[b] *= p;

			if (w_compute_bound > 0)
			{
				v = theta[0];
				for (b=1; b<K; b++) if (v > theta[b]) v = theta[b];
				LB += v * w_compute_bound;
			}

			for (e=B->first_in; e; e=e->next_in)
			{
				if (e->weight_forward == 1)
				{
					for (b=0; b<K; b++) e->m[b] -= theta[b];
				}
				else
				{
					for (b=0; b<K; b++) e->m[b] -= e->weight_forward*theta[b];
				}
			}
		}

		if (compute_solution)
		{
			cost_fw = ComputeCost();
		}

		if (options.verbose)
		{
			if (options.print_times) printf("iter %d [%.3f secs]: ", iter+1, get_time() - time_start);
			else                     printf("iter %d: ", iter+1);
			if (compute_solution) printf("lower bound=%f, cost_fw=%f\n", LB, cost_fw);
			else                  printf("lower bound=%f\n", LB);
		}

#ifdef SRMP_VNK_DEBUG
		double LB_check = ComputeLowerBound();
		assert (fabs(LB - LB_check) < 1e-5);
#endif
	}

	delete [] theta;
	delete [] restricted_buf;
	delete [] seq.arr;

	return LB;
}



