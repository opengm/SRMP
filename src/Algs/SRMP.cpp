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
#include <limits>
#include <srmp/Algs/util.h>
#include <srmp/SRMP.h>

using namespace srmpLib;


double Energy::InitSRMP(Sequence& seq, Options& options)
{
	int t;
	Node* i;
	Factor* A;
	NonSingletonFactor* _A;
	Edge* e;

	// set 'separator_num' and 'order'
	int separator_num = 0;
	double LB_init = 0;
	for (i=nodes; i<nodes+node_num; i++) { if (!i->is_removed) { i->tmp1 = i->tmp2 = 0; separator_num ++; } }
	for (_A=factors.ScanFirst(); _A; _A=factors.ScanNext())
	{
		if (_A->is_removed) continue;
		_A->tmp1 = _A->tmp2 = 0;
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
	// use Factor::rep for storing order (temporarily)
	seq.arr[0].A->rep = NULL;
	for (t=1; t<separator_num; t++) seq.arr[t].A->rep = seq.arr[t-1].A->rep+1;

	seq.arity_max = 0;
	seq.K_max = 0;
	for (t=0; t<separator_num; t++)
	{
		A = seq.arr[t].A;

		if (seq.arity_max < A->arity) seq.arity_max = A->arity;
		if (seq.K_max < A->K) seq.K_max = A->K;

		if (A->tmp1)
		{
			A->compute_bound = (A->arity == 1) ? 1 : 0;
		}
		else
		{
			A->tmp1 = 1;
			A->compute_bound = 1;
		}
		for (e=A->first_in; e; e=e->next_in)
		{
			if (e->A->tmp1)
			{
				e->is_bw = 1;
				e->compute_bound = 0;
			}
			else
			{
				e->is_bw = 0; 
				e->A->tmp1 = 1;
				e->compute_bound = 1;
			}
		}
	}
	for (t=separator_num-1; t>=0; t--)
	{
		A = seq.arr[t].A;
		A->tmp2 = 1;
		for (e=A->first_in; e; e=e->next_in)
		{
			if (e->A->tmp2) e->is_fw = 1;
			else { e->is_fw = 0; e->A->tmp2 = 1; }
		}
	}

	// set weights
	int default_weight = 1;
	if (options.TRWS_weighting != (double)((int)options.TRWS_weighting)) default_weight = 10;
	for (t=0; t<separator_num; t++)
	{
		A = seq.arr[t].A;

		////////////////////////////////////////////////////
		int w_forward_out = 0, w_backward_out = 0;
		if (A->arity > 1)
		{
			for (e=((NonSingletonFactor*)A)->first_out; e; e=e->next_out)
			{
				if (e->B->rep > A->rep) w_forward_out += default_weight;
				else                    w_backward_out += default_weight;
			}
		}

		////////////////////////////////////////////////////
		int w, w_forward_in = 0, w_backward_in = 0, w_total_in = 0;
		

		for (e=A->first_in; e; e=e->next_in)
		{
			int e_weight = default_weight; // can be changed to other values (perhaps, dependent on arities?)

			if (e->is_fw)
			{
				e->weight_forward = e_weight;
				w_forward_in += e_weight;
			}
			else e->weight_forward = 0;

			if (e->is_bw)
			{
				e->weight_backward = e_weight;
				w_backward_in += e_weight;
			}
			else e->weight_backward = 0;

			w_total_in += e_weight;
		}

		// if TRWS_weighting == 0 then A->weight_forward = w_forward_out + w_forward_in
		// if TRWS_weighting == 1 then A->weight_forward = w_forward_out + max ( w_forward_in, w_total_in-w_forward_in )

		int delta_forward = (w_total_in - w_forward_in) - w_forward_in; if (delta_forward < 0) delta_forward = 0;
		//int delta_forward = w_backward_in - w_forward_in; if (delta_forward < 0) delta_forward = 0;
		delta_forward = (int)(options.TRWS_weighting*delta_forward);
		w = w_forward_out + w_forward_in + delta_forward;
		A->weight_forward = (unsigned)w;
		if ((int)A->weight_forward != w) { printf("Error: capacity of Factor::weight_forward is not enough!\n"); exit(1); }

		// similarly for backward

		int delta_backward = (w_total_in - w_backward_in) - w_backward_in; if (delta_backward < 0) delta_backward = 0;
		//int delta_backward = w_forward_in - w_backward_in; if (delta_backward < 0) delta_backward = 0;
		delta_backward = (int)(options.TRWS_weighting*delta_backward);
		w = w_backward_out + delta_backward + w_backward_in;
		A->weight_backward = (unsigned)w;
		if ((int)A->weight_backward != w) { printf("Error: capacity of Factor::weight_backward is not enough!\n"); exit(1); }

		/////////

		if (A->weight_forward + w_forward_in == 0) A->weight_forward = 1;
		if (A->weight_backward + w_backward_in == 0) A->weight_backward = 1;
	}

	return LB_init;
}



double Energy::SolveSRMP(Options& options)
{
	double time_start = get_time();

	Sequence seq;
	double LB_init = InitSRMP(seq, options);

	Node* i;
	int t, iter, iter_max = options.iter_max;
	Factor* B;
	Edge* e;

	double v, LB = 0, LB_prev = 0;
	double cost_fw, cost_bw;

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
		for (t=0; t<seq.num; t++)
		{
			B = seq.arr[t].A;
			int b, K = B->K;
			if (B->arity == 1) memcpy(theta, B->data, K*sizeof(double));
			else SRMP_COMPUTE_PARTIAL_REPARAMETERIZATION((NonSingletonFactor*)B, theta);

			int w_total = B->weight_forward;
			for (e=B->first_in; e; e=e->next_in)
			{
				if (e->is_bw)
				{
				    SRMP_SEND_MESSAGE(e);
				}
				for (b=0; b<K; b++) theta[b] += e->m[b];
			}

			if (compute_solution)
			{
				ComputeSolution(B, restricted_buf);
			}

			double p = 1.0 / w_total;
			for (b=0; b<K; b++) theta[b] *= p;

			for (e=B->first_in; e; e=e->next_in)
			{
				if (e->weight_forward == 0) continue;
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

		// backward pass
		if (compute_solution)
		{
			for (i=nodes; i<nodes+node_num; i++) i->solution = -1;
		}
		LB = LB_init;
		for (t=seq.num-1; t>=0; t--)
		{
			B = seq.arr[t].A;
			int b, K = B->K;

			if (B->arity == 1) memcpy(theta, B->data, K*sizeof(double));
			else SRMP_COMPUTE_PARTIAL_REPARAMETERIZATION((NonSingletonFactor*)B, theta);

			int B_weight = B->weight_backward;
			for (e=B->first_in; e; e=e->next_in)
			{
				if (e->is_fw || e->compute_bound)
				{
					v = SRMP_SEND_MESSAGE(e);
					if (e->compute_bound) LB += v;
				}
				for (b=0; b<K; b++) theta[b] += e->m[b]; 
				B_weight -= e->weight_backward;
			}

			if (compute_solution)
			{
				ComputeSolution(B, restricted_buf);
			}

			double p = 1.0 / B->weight_backward;
			for (b=0; b<K; b++) theta[b] *= p;

			if (B->compute_bound && B->weight_backward > 0)
			{
				v = theta[0];
				for (b=1; b<K; b++) if (v > theta[b]) v = theta[b];
				LB += v * B_weight;
			}

			for (e=B->first_in; e; e=e->next_in)
			{
				if (e->weight_backward == 0) continue;
				if (e->weight_backward == 1)
				{
					for (b=0; b<K; b++) e->m[b] -= theta[b];
				}
				else
				{
					for (b=0; b<K; b++) e->m[b] -= e->weight_backward*theta[b];
				}
			}
		}
		if (compute_solution)
		{
			cost_bw = ComputeCost();
		}

		if (options.verbose)
		{
			if (options.print_times) printf("iter %d [%.3f secs]: ", iter+1, get_time() - time_start);
			else                     printf("iter %d: ", iter+1);
			if (compute_solution) printf("lower bound=%f, cost_fw=%f, cost_bw=%f\n", LB, cost_fw, cost_bw);
			else                 printf("lower bound=%f\n", LB);
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







double Energy::Solve(Options& options)
{
	if (dual_graph) return SolveDualGraph(options);

	if (options.verbose)
	{
		if (options.method == Options::SRMP)
		{
			printf("SRMP (TRW-S weighting=%f)\n", options.TRWS_weighting);
		}
		else if (options.method == Options::MPLP)
		{
			printf("MPLP\n");
		}
		else if (options.method == Options::MPLP_BW)
		{
			printf("MPLP_BW\n");
		}
		else if (options.method == Options::CMP)
		{
			printf("CMP\n");
		}
		else
		{
			printf("Unknown method\n");
			exit(1);
		}
		printf("  iter_max=%d\n", options.iter_max);
		printf("  time_max=%.3f secs\n", options.time_max);
		printf("  eps=%.15f\n", options.eps);
		printf("  compute solution every %d iters\n", options.compute_solution_period);
		printf("  sort_flag=%d\n", options.sort_flag);
	}

	double time_start = get_time();

	if (!edges) SetFullEdges();
	InitEdges();

	double LB;

	if (options.iter_max == 0)
	{
		if (options.verbose)
		{
			LB = ComputeLowerBound();
			printf("lower_bound=%f\n", LB);
			return LB;
		}
		return -std::numeric_limits<double>::infinity();
	}

	if      (options.method == Options::SRMP) LB = SolveSRMP(options);
	else if (options.method == Options::CMP)  LB = SolveCMP(options);
	else                                      LB = SolveMPLP(options);

	if (options.verbose) printf("Done (%.3f secs). cost_best=%f\n", get_time()-time_start, cost_best);

	return LB;
}

