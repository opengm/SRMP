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
#include "../FactorTypes/PairwiseType.h"
#include "../FactorTypes/GeneralType.h"


Energy::Energy(int _node_num_max)
	: node_num(0), node_num_max(_node_num_max), arity_max(1),
	  factors(512), edges(NULL), buf(4096), is_solution_best_initialized(false), is_cost_best_valid(false),
	  factor_type_pairwise(NULL), factor_type_general(NULL), 
	  primal_graph(NULL), factor_type_pairwise_dual(NULL), dual_graph(NULL), dual_sequence(NULL)
{
	nodes = new Node[node_num_max];
}

Energy::~Energy()
{
	delete [] nodes;
	if (edges) delete edges;
	if (factor_type_pairwise) delete factor_type_pairwise;
	if (factor_type_general) delete factor_type_general;

	if (factor_type_pairwise_dual) delete factor_type_pairwise_dual;
	if (dual_graph) delete dual_graph;
	if (dual_sequence) delete dual_sequence;
}

Energy::NodeId Energy::AddNode(int K, double* costs)
{
	Node* i = nodes + node_num;
	i->arity = 1;
	i->K = K;
	if (costs) i->data = costs;
	else
	{
		i->data = (double*) buf.Alloc(K*sizeof(double));
		memset(i->data, 0, K*sizeof(double));
	}
	i->first_in = NULL;
	i->is_removed = 0;
	return node_num ++;
}

Energy::FactorId Energy::AddUnaryFactor(NodeId _i, double* costs)
{
	int k;
	Node* i = nodes + _i;
	if (costs)
	{
		for (k=0; k<i->K; k++) i->data[k] += costs[k];
	}
	return i;
}

Energy::FactorId Energy::AddPairwiseFactor(NodeId _i, NodeId _j, double* costs)
{
	NodeId node_indexes[2] = { _i, _j };
	return AddFactor(2, node_indexes, costs);
}

Energy::FactorId Energy::AddFactor(int arity, NodeId* node_indexes, double* costs, FactorType* type, unsigned flags)
{
	if (arity == 1) return AddUnaryFactor(node_indexes[0], costs);

	int i;
	NonSingletonFactor* A = factors.New();
	A->arity = arity;
	A->is_unsorted = 0;
	A->is_removed = 0;
	for (i=0; i<arity-1; i++)
	{
		if (node_indexes[i] > node_indexes[i+1]) { A->is_unsorted = 1; break; }
	}
	A->K = 1;
	A->nodes = (Node**) buf.Alloc((1+A->is_unsorted)*arity*sizeof(Node*));
	for (i=0; i<arity; i++)
	{
		A->nodes[i] = nodes + node_indexes[i];
		A->K *= A->nodes[i]->K;
	}
	if (A->is_unsorted)
	{
		for (i=0; i<arity; i++) A->nodes[i+arity] = A->nodes[i];
		quickSort<Node*>(A->nodes+arity, 0, arity-1);
	}
	A->first_in = A->first_out = NULL;

	if (arity_max < arity) arity_max = arity;

	InitFactor(A, costs, type);

	return A;
}




void Energy::InitFactor(NonSingletonFactor* A, double* user_data, FactorType* type, unsigned flags)
{
	if (type) A->type = type;
	else
	{
		if (A->arity == 2)
		{
			if (!factor_type_pairwise)
			{
				factor_type_pairwise = new PairwiseFactorType;
			}
			A->type = factor_type_pairwise;
		}
		else
		{
			if (!factor_type_general)
			{
				factor_type_general = new GeneralFactorType;
			}
			A->type = factor_type_general;
		}
	}
	A->type->InitFactor(A, user_data, flags);
}


void Energy::GetFactorCosts(NonSingletonFactor* A, double* costs, void* _buf)
{
	int i, n = A->arity;
	Node** Anodes = A->nodes;

	int* solution_old = (int*)_buf;

	int k = 0;
	for (i=0; i<n; i++)
	{
		solution_old[i] = Anodes[i]->solution;
		Anodes[i]->solution = 0;
	}
	costs[k ++] = A->type->GetCost(A);
	while ( 1 )
	{
		for (i=n-1; i>=0; i--)
		{
			if (Anodes[i]->solution < Anodes[i]->K-1) break;
			Anodes[i]->solution = 0;
		}
		if (i<0) break;
		Anodes[i]->solution ++;
		costs[k ++] = A->type->GetCost(A);
	}
	for (i=0; i<n; i++)
	{
		Anodes[i]->solution = solution_old[i];
	}
}

void Energy::ConvertFactorToGeneral(NonSingletonFactor* A, void* _buf)
{
	double* costs = (double*) _buf;
	GetFactorCosts(A, costs, (void*)(costs + A->K));

	InitFactor(A, costs, NULL);

	Edge* e;
	for (e=A->first_out; e; e=e->next_out)
	{
		if (e->m)
		{
			A->type->InitEdge(e);
		}
	}
}


void Energy::InitEdges()
{
	NonSingletonFactor* A;
	Edge* e;
	ReusableBuffer convert_factor_rbuf;

	int messages_size = 0;
	for (e=edges->ScanFirst(); e; e=edges->ScanNext())
	{
		if (!e->m)
		{
			messages_size += e->B->K;
			e->A->type->InitEdge(e);
		}
	}

	double* messages = (double*) buf.Alloc(messages_size*sizeof(double));
	memset(messages, 0, messages_size*sizeof(double));

	for (e=edges->ScanFirst(); e; e=edges->ScanNext())
	{
		if (!e->m)
		{
			e->m = messages;
			messages += e->B->K;
		}
	}

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (!A->type->PrepareFactor(A))
		{
			void* _buf = convert_factor_rbuf.Alloc(A->K*sizeof(double) + A->arity*sizeof(int));
			ConvertFactorToGeneral(A, _buf);
		}
	}

	if (is_solution_best_initialized) // SolveXXX() is called second time. New factors may have been added, so recompute cost of best solution
	{
		is_cost_best_valid = false;
		Node* i;
		for (i=nodes; i<nodes+node_num; i++) i->solution = i->solution_best;
		double cost = ComputeCost();
		double LB = ComputeLowerBound();
		printf("start: lower bound=%f, cost=%f\n", LB, cost);
	}
}

double Energy::ComputeCost()
{
	double cost = 0;
	Node* i;
	NonSingletonFactor* A;
	for (i=nodes; i<nodes+node_num; i++)
	{
		cost += i->data[i->solution];
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		cost += A->type->GetCost(A);
	}
	if (!is_cost_best_valid || cost_best > cost)
	{
		is_cost_best_valid = true;
		is_solution_best_initialized = true;
		cost_best = cost;
		for (i=nodes; i<nodes+node_num; i++) i->solution_best = i->solution;
	}
	if (primal_graph && cost >= ENERGY_INFTY/2) // we are in the dual energy; the solution is inconsistent!
	{
		cost = primal_graph->ConvertSolutionDualToPrimal();
		primal_graph->dual_solution_was_inconsistent = true;
	}
	return cost;
}

void Energy::AddTriplet(NodeId _i, NodeId _j, NodeId _k, bool add_all_edges)
{
	if (_i<0 || _j<0 || _k<0 || _i>=node_num || _j>=node_num || _k>=node_num || _i==_j || _i==_k || _j==_k)
	{
		printf("Incorrect arguments to AddTriplet()\n");
		exit(1);
	}

	int r;
	FactorId f[3]; // fij, fjk, fik;
	NonSingletonFactor* A;
	Edge* e;

	for (r=0; r<3; r++)
	{
		Node* i;
		Node* j;
		if (r==0)      { i = nodes+_i; j = nodes+_j; }
		else if (r==1) { i = nodes+_j; j = nodes+_k; }
		else           { i = nodes+_i; j = nodes+_k; }

		for (e=i->first_in; e; e=e->next_in)
		{
			A = e->A;
			if (A->arity == 2 && (A->nodes[0] == j || A->nodes[1] == j))
			{
				f[r] = A;
				break;
			}
		}
		if (!e)
		{
			f[r] = AddPairwiseFactor((int)(i-nodes), (int)(j-nodes), NULL);
			if (add_all_edges)
			{
				AddRelaxationEdge(f[r], i);
				AddRelaxationEdge(f[r], j);
			}
		}
	}
	NodeId arr[3] = { _i, _j, _k };
	FactorId ijk = AddFactor(3, arr, NULL);
	AddRelaxationEdge(ijk, f[0]);
	AddRelaxationEdge(ijk, f[1]);
	AddRelaxationEdge(ijk, f[2]);
}




void Energy::ComputeSolution(Factor* B, void* _buf)
{
	int i, b, K = B->K;
	Edge* e;
	double* theta = (double*)_buf;
	double* m = theta + B->K;
	int* ComputeRestrictedMinimum_buf = (int*)(m + B->K); // of size 4*B->arity

	bool exists_unlabeled;
	if (B->arity == 1) exists_unlabeled = (((Node*)B)->solution < 0);
	else
	{
		for (i=0; i<B->arity; i++) if (((NonSingletonFactor*)B)->nodes[i]->solution < 0) break;
		exists_unlabeled = (i < B->arity);
	}
	if (exists_unlabeled)
	{
		if (B->arity == 1) memcpy(theta, B->data, K*sizeof(double));
		else COMPUTE_PARTIAL_REPARAMETERIZATION((NonSingletonFactor*)B, theta);

		for (e=B->first_in; e; e=e->next_in)
		{
			int unlabeled_num = 0;
			for (i=0; i<e->A->arity; i++) if (e->A->nodes[i]->solution >= 0) unlabeled_num ++;
			if (unlabeled_num > 0 && unlabeled_num < e->A->arity)
			{
				double* m_old = e->m;
				e->m = m;
				SEND_RESTRICTED_MESSAGE(e);
				e->m = m_old;
				for (b=0; b<K; b++) theta[b] += m[b];
			}
			else
			{
				for (b=0; b<K; b++) theta[b] += e->m[b];
			}
		}
		if (B->arity == 1)
		{
			int b_best = 0;
			double v = theta[0];
			for (b=1; b<K; b++) if (v > theta[b]) { b_best = b; v = theta[b]; }
			((Node*)B)->solution = b_best;
		}
		else
		{
			((NonSingletonFactor*)B)->ComputeRestrictedMinimum(theta, ComputeRestrictedMinimum_buf);
		}
	}
}


void Energy::NonSingletonFactor::ComputeRestrictedMinimum(double* theta, int* _buf)
{
	int i, k = 0, n = 0;
	int* Kfactor_array = _buf;
	int* K_array = _buf + arity;
	int* index_array = _buf + 2*arity;
	int* labeling = _buf + 3*arity;

	int K_factor = 1;
	for (i=arity-1; i>=0; i--)
	{
		if (nodes[i]->solution >= 0) k += nodes[i]->solution * K_factor;
		else
		{
			nodes[i]->solution = 0;
			K_array[n] = nodes[i]->K;
			Kfactor_array[n] = K_factor;
			index_array[n] = i;
			labeling[n ++] = 0;
		}
		K_factor *= nodes[i]->K;
	}

	if (n == arity)
	{
		// everything is unlabeled
		int k_best = 0;
		double v_best = theta[0];
		for (k=1; k<K; k++)
		{
			if (v_best > theta[k]) { k_best = k; v_best = theta[k]; }
		}
		for (i=arity-1; ;i--)
		{
			nodes[i]->solution = k_best % nodes[i]->K;
			if (i == 0) return;
			k_best /= nodes[i]->K;
		}
	}

	double v_best = theta[k];
	while ( 1 )
	{
		for (i=0; i<n; i++)
		{
			if (labeling[i] < K_array[i]-1) break;
			k -= labeling[i]*Kfactor_array[i];
			labeling[i] = 0;
		}
		if (i==n) break;
		labeling[i] ++;
		k += Kfactor_array[i];
		if (v_best > theta[k])
		{
			v_best = theta[k];
			for (i=0; i<n; i++)
			{
				nodes[index_array[i]]->solution = labeling[i];
			}
		}
	}
}

void Energy::NonSingletonFactor::ComputeRestrictedMinimum(double* theta, Factor* B, double* thetaB, int* _buf)
{
	int i, k = 0, n = 0;
	int* Kfactor_array = _buf;
	int* K_array = _buf + arity;
	int* labeling = _buf + 2*arity;

	int kB = 0, nB = 0;
	int* KBfactor_array = _buf + 3*arity;
	Node** Bnodes = (B->arity == 1) ? (Node**)(&B) : ((NonSingletonFactor*)B)->nodes;

	int K_factor = 1;
	for (i=arity-1; i>=0; i--)
	{
		int iB, KB_factor = 1;
		for (iB=B->arity-1; iB>=0; iB--)
		{
			if (nodes[i] == Bnodes[iB]) break;
			KB_factor *= Bnodes[iB]->K;
		}
		if (iB < 0) KB_factor = 0;

		if (nodes[i]->solution >= 0)
		{
			kB += nodes[i]->solution * KB_factor;
			k += nodes[i]->solution * K_factor;
		}
		else
		{
			KBfactor_array[n] = KB_factor;
			K_array[n] = nodes[i]->K;
			Kfactor_array[n] = K_factor;
			labeling[n ++] = 0;
		}
		K_factor *= nodes[i]->K;
	}

	for (i=0; i<B->K; i++) thetaB[i] = 1e100;
	thetaB[kB] = theta[k];
	while ( 1 )
	{
		for (i=0; i<n; i++)
		{
			if (labeling[i] < K_array[i]-1) break;
			k -= labeling[i]*Kfactor_array[i];
			kB -= labeling[i]*KBfactor_array[i];
			labeling[i] = 0;
		}
		if (i==n) break;
		labeling[i] ++;
		k += Kfactor_array[i];
		kB += KBfactor_array[i];
		if (thetaB[kB] > theta[k])
		{
			thetaB[kB] = theta[k];
		}
	}
}


double Energy::ComputeLowerBound()
{
	Node* i;
	NonSingletonFactor* A;
	Edge* e;
	int b;
	double LB = 0;

	ReusableBuffer theta_rbuf;

	for (i=nodes; i<nodes+node_num; i++)
	{
		double* theta = (double*) theta_rbuf.Alloc(i->K*sizeof(double));
		memcpy(theta, i->data, i->K*sizeof(double));
		for (e=i->first_in; e; e=e->next_in)
		{
			for (b=0; b<i->K; b++) theta[b] += e->m[b];
		}
		double v_min = theta[0];
		for (b=1; b<i->K; b++)
		{
			if (v_min > theta[b]) v_min = theta[b];
		}
		LB += v_min;
	}

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		if (!A->first_out)
		{
			double* Arep_old = A->rep;
			A->rep = NULL;
			LB += SEND_MPLP_MESSAGES(A, false);
			A->rep = Arep_old;
			continue;
		}

		e = A->first_out;
		double* m_old = e->m;
		double* m_new = e->m = (double*) theta_rbuf.Alloc(e->B->K*sizeof(double));
		memcpy(m_new, m_old, e->B->K*sizeof(double));
		LB += SEND_MESSAGE(e);
		e->m = m_old;
		double v_min = -m_old[0] + m_new[0];
		for (b=1; b<e->B->K; b++)
		{
			if (v_min > -m_old[b] + m_new[b]) v_min = -m_old[b] + m_new[b];
		}
		LB += v_min;

		A->type->PrepareFactor(A);
	}

	return LB;
}

void Energy::PrintStats()
{
	NonSingletonFactor* A;
	Edge* e;
	int k;
	int* stat = new int[arity_max + arity_max*arity_max];
	for (k=0; k<arity_max + arity_max*arity_max; k++) stat[k] = 0;

	stat[0] = node_num;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		stat[A->arity-1] ++;
		for (e=A->first_out; e; e=e->next_out) stat[arity_max + (A->arity-1)*arity_max + e->B->arity-1] ++;
	}

	printf("Factor stats:\n");
	for (k=0; k<arity_max; k++)
	{
		printf("  %d: %d\n", k+1, stat[k]);
	}
	printf("Edge stats:\n");
	for (k=0; k<arity_max*arity_max; k++)
	{
		if (stat[arity_max+k] == 0) continue;
		printf("  %d->%d: %d\n", (k/arity_max)+1, (k%arity_max)+1, stat[arity_max+k]);
	}

	delete [] stat;
}


#ifdef _MSC_VER
#pragma warning(disable: 4996) /* Disable deprecation */
#endif

void Energy::SaveUAI(char* filename, bool sort_factors, bool save_reparameterization)
{
	Node* i;
	Factor* A;
	ReusableBuffer rbuf;

	int t, k, factor_num = node_num;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		factor_num ++;
	}

	FactorPtr* F = new FactorPtr[factor_num];

	factor_num = 0;
	for (i=nodes; i<nodes+node_num; i++) F[factor_num ++].A = i;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		F[factor_num ++].A = A;
	}

	if (sort_factors)
	{
		quickSort<FactorPtr>(F, 0, factor_num-1);
	}

	///////////////////////////////////////////////////////////

	FILE* fp = fopen(filename, "w");
	if (!fp) { printf("Can't open %s for writing\n", filename); exit(1); }

	fprintf(fp, "MARKOV\n");
	fprintf(fp, "%d\n", node_num);
	for (i=nodes; i<nodes+node_num; i++)
	{
		fprintf(fp, "%d ", i->K);
	}
	fprintf(fp, "\n");
	fprintf(fp, "%d\n", factor_num);

	for (t=0; t<factor_num; t++)
	{
		A = F[t].A;
		fprintf(fp, "%d", A->arity);
		if (A->arity == 1) fprintf(fp, " %d", (int)(((Node*)A) - nodes));
		else
		{
			for (k=0; k<A->arity; k++) fprintf(fp, " %d", (int)( ((NonSingletonFactor*)A)->nodes[k] - nodes ));
		}
		fprintf(fp, "\n");
	}
	if (!save_reparameterization)
	{
		for (t=0; t<factor_num; t++)
		{
			A = F[t].A;
			fprintf(fp, "\n%d\n", A->K);
			double* costs;
			if (A->arity == 1) costs = A->data;
			else
			{
				costs = (double*) rbuf.Alloc(A->K*sizeof(double) + A->arity*sizeof(int));
				void* _buf = costs + A->K;
				GetFactorCosts((NonSingletonFactor*)A, costs, _buf);
			}
			for (k=0; k<A->K; k++) fprintf(fp, "%f ", -costs[k]);
			fprintf(fp, "\n");
		}
	}
	else
	{
		for (t=0; t<factor_num; t++)
		{
			A = F[t].A;
			fprintf(fp, "\n%d\n", A->K);

			double* costs = (double*) rbuf.Alloc(A->K*sizeof(double));

			if (A->arity == 1) memcpy(costs, A->data, A->K*sizeof(double));
			else
			{
				COMPUTE_PARTIAL_REPARAMETERIZATION((NonSingletonFactor*)A, costs);
			}

			Edge* e;
			for (e=A->first_in; e; e=e->next_in)
			{
				for (k=0; k<A->K; k++) costs[k] += e->m[k];
			}

			for (k=0; k<A->K; k++) fprintf(fp, "%f ", -costs[k]);
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
	///////////////////////////////////////////////////////////

	delete [] F;
}
