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

double Energy::InitMPLP(Energy::Sequence& seq, Options& options)
{
	int t;
	Node* i;
	NonSingletonFactor* A;
	Edge* e;

	seq.num = 0;
	double LB_init = 0;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		if (A->first_out) seq.num ++;
		else if (!A->first_in)
		{
			A->rep = NULL;
			LB_init += SEND_MPLP_MESSAGES(A, false);
		}
	}
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (i->is_removed) continue;
		if (!i->first_in)
		{
			int k;
			double v = i->data[0];
			i->solution = 0;
			for (k=1; k<i->K; k++)
			{
				if (v > i->data[k]) { v = i->data[k]; i->solution = k; }
			}
			LB_init += v;
		}
	}

	seq.arr = new FactorPtr[seq.num];

	t = 0;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (!A->is_removed && A->first_out) seq.arr[t ++].A = A;
	}

	SortSequence(seq, options);

	if (options.method == Options::MPLP_BW)
	{
		for (t=seq.num-1; t>=0; t--)
		{
			A = (NonSingletonFactor*)seq.arr[t].A;
			A->compute_bound = 1;
			for (e=A->first_out; e; e=e->next_out) e->B->compute_bound = 0;
		}
	}
	else
	{
		for (t=0; t<seq.num; t++)
		{
			A = (NonSingletonFactor*)seq.arr[t].A;
			A->compute_bound = 1;
			for (e=A->first_out; e; e=e->next_out) e->B->compute_bound = 0;
		}
	}

	return LB_init;
}
	
	
double Energy::SolveMPLP(Options& options)
{
	double time_start = get_time();

	Sequence seq;
	double LB_init = InitMPLP(seq, options);

	Node* i;
	int b, t, iter, iter_max = options.iter_max;
	NonSingletonFactor* A;
	Edge* e;

	double v, LB = 0, LB_prev = 0;
	double cost_fw = 0, cost_bw = 0;

	//////////////////////////////////////////

	int rep_size = 0;
	for (i=nodes; i<nodes+node_num; i++)
	{
		rep_size += i->K;
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->first_in) rep_size += A->K;
	}

	double* rep_array = new double[rep_size];
	rep_size = 0;
	for (i=nodes; i<nodes+node_num; i++)
	{
		i->rep = rep_array + rep_size; rep_size += i->K;
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->first_in) { A->rep = rep_array + rep_size; rep_size += A->K; }
		else A->rep = NULL;

		A->weight_forward = 0;
		for (e=A->first_out; e; e=e->next_out) e->weight_forward = 1;
	}

	//////////////////////////////////////////

	for (iter=0; iter<iter_max; iter++)
	{
		if ((iter > 1 && LB<LB_prev + options.eps) || get_time()-time_start > options.time_max)
			iter_max = iter+1; // last iteration, with compute_solution=true
		LB_prev = LB;

		bool compute_solution = false;
		if (iter != 0 && iter % options.compute_solution_period == 0) compute_solution = true;
		if (iter == iter_max-1) compute_solution = true;

		if ( iter % 10 == 0 )
		{
			// recompute Factor::rep
			for (i=nodes; i<nodes+node_num; i++)
			{
				memcpy(i->rep, i->data, i->K*sizeof(double));
				for (e=i->first_in; e; e=e->next_in)
				{
					for (b=0; b<i->K; b++) i->rep[b] += e->m[b];
				}
			}
			for (A=factors.ScanFirst(); A; A=factors.ScanNext())
			{
				if (!A->rep) continue;
				A->type->ComputePartialReparameterization(A, A->rep);
				for (e=A->first_in; e; e=e->next_in)
				{
					for (b=0; b<A->K; b++) A->rep[b] += e->m[b];
				}
			}
		}

		// forward pass
		LB = LB_init;
		if (compute_solution)
		{
			for (i=nodes; i<nodes+node_num; i++)
			{
				if (i->first_in) i->solution = -1;
			}
		}
		for (t=0; t<seq.num; t++)
		{
			A = (NonSingletonFactor*) seq.arr[t].A;

			for (e=A->first_out; e; e=e->next_out)
			{
				for (b=0; b<e->B->K; b++)
				{
					e->B->rep[b] -= e->m[b];
					e->m[b] = -e->B->rep[b];
				}
			}

			v = SEND_MPLP_MESSAGES(A, compute_solution);
			if (A->compute_bound) LB += v;

			for (e=A->first_out; e; e=e->next_out)
			{
				for (b=0; b<e->B->K; b++)
				{
					e->m[b] += e->B->rep[b];
				}
			}
		}
		if (compute_solution)
		{
			cost_fw = ComputeCost();
		}

		// backward pass
		if (options.method == Options::MPLP_BW)
		{
			LB = LB_init;
			if (compute_solution)
			{
				for (i=nodes; i<nodes+node_num; i++)
				{
					if (i->first_in) i->solution = -1;
				}
			}
			for (t=seq.num-1; t>=0; t--)
			{
				A = (NonSingletonFactor*) seq.arr[t].A;

				for (e=A->first_out; e; e=e->next_out)
				{
					for (b=0; b<e->B->K; b++)
					{
						e->B->rep[b] -= e->m[b];
						e->m[b] = -e->B->rep[b];
					}
				}

				v = SEND_MPLP_MESSAGES(A, compute_solution);
				if (A->compute_bound) LB += v;

				for (e=A->first_out; e; e=e->next_out)
				{
					for (b=0; b<e->B->K; b++)
					{
						e->m[b] += e->B->rep[b];
					}
				}
			}
			if (compute_solution)
			{
				cost_bw = ComputeCost();
			}
		}

		if (options.verbose)
		{
			if (options.print_times) printf("iter %d [%.3f secs]: ", iter+1, get_time() - time_start);
			else                     printf("iter %d: ", iter+1);
			if (compute_solution)
			{
				if (options.method == Options::MPLP_BW) printf("lower bound=%f, cost_fw=%f, cost_bw=%f\n", LB, cost_fw, cost_bw);
				else                                    printf("lower bound=%f, cost_fw=%f\n", LB, cost_fw);
			}
			else                                        printf("lower bound=%f\n", LB);
		}

#ifdef VNK_DEBUG
		double LB_check = ComputeLowerBound();
		assert (fabs(LB - LB_check) < 1e-5);
#endif
	}

	delete [] rep_array;

	return LB;
}


