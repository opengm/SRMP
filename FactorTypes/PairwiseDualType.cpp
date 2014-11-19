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
#include "PairwiseDualType.h"
#include "GeneralType.h"
#include "../Algs/util.h"

#define SRMP_WHICH_NODE user1

int PairwiseDualFactorType::ComputeRestriction(Energy::Edge* e, int k)
{
	int i;
	for (i=e->A->arity-1; i>=0; i--)
	{
		int K = e->A->nodes[i]->K;
		e->A->nodes[i]->solution = k % K;
		k /= K;
	}
	if (e->B->arity == 1) return ((Energy::Node*)e->B)->solution;

	k = 0;
	int Kfactor = 1;
	for (i=e->B->arity-1; i>=0; i--)
	{
		int K = ((Energy::NonSingletonFactor*)e->B)->nodes[i]->K;
		k += Kfactor * ((Energy::NonSingletonFactor*)e->B)->nodes[i]->solution;
		Kfactor *= K;
	}
	return k;
}

PairwiseDualFactorType::PairwiseDualFactorType() : buf(4096)
{
}

PairwiseDualFactorType::~PairwiseDualFactorType()
{
}

void PairwiseDualFactorType::InitFactor(Energy::NonSingletonFactor* A, double* user_data, unsigned flags)
{
	Energy::Edge* e[2] = { ((Energy::Edge**)user_data)[0], ((Energy::Edge**)user_data)[1] };
	assert(A->arity == 2);
	assert(e[0]->B == e[1]->B && e[0]->A->K == A->nodes[0]->K && e[1]->A->K == A->nodes[1]->K);

	A->data = (double*) buf.Alloc(2*sizeof(Energy::Edge*));
	memcpy(A->data, user_data, 2*sizeof(Energy::Edge*));
}

double PairwiseDualFactorType::GetCost(Energy::NonSingletonFactor* AB)
{
	Energy::Edge* e[2] = { ((Energy::Edge**)AB->data)[0], ((Energy::Edge**)AB->data)[1] };
	Energy::NonSingletonFactor* A[2] = { e[0]->A, e[1]->A };
	int i, k;

	for (i=0; i<A[1]->arity; i++) A[1]->nodes[i]->solution = -1;

	k = AB->nodes[0]->solution;
	for (i=A[0]->arity-1; i>=0; i--)
	{
		int Ki = A[0]->nodes[i]->K;
		A[0]->nodes[i]->solution = k % Ki;
		k /= Ki;
	}

	k = AB->nodes[1]->solution;
	for (i=A[1]->arity-1; i>=0; i--)
	{
		int Ki = A[1]->nodes[i]->K;
		if (A[1]->nodes[i]->solution >= 0 && A[1]->nodes[i]->solution != k % Ki) return SRMP_ENERGY_INFTY;
		k /= Ki;
	}

	return 0;
}

void PairwiseDualFactorType::InitEdge(Energy::Edge* e)
{
	Energy::NonSingletonFactor* A = e->A;

	if ((Energy::Node*)e->B == A->nodes[0]) e->SRMP_WHICH_NODE = 0;
	else                                    e->SRMP_WHICH_NODE = 1;
}

bool PairwiseDualFactorType::PrepareFactor(Energy::NonSingletonFactor* A)
{
	if (A->first_in) return false; // if there are incoming edges then convert factor to general 'PairwiseType'
	if (!A->first_out || !A->first_out->next_out || A->first_out->next_out->next_out
		|| A->first_out->B == A->first_out->next_out->B) return false; // degenerate cases, probably will never occur
	return true;
}

void PairwiseDualFactorType::ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* theta)
{
	// this function should never be called since incoming edges are not allowed (as specified in PrepareFactor())
	printf("Error: ComputePartialReparameterization() should not be called for this type. (Trying to save non-standard factor?");
	exit(1);
}


double PairwiseDualFactorType::SendMessage(Energy::Edge* e0)
{
	Energy::NonSingletonFactor* AC = e0->A;
	Energy::Edge* eA = ((Energy::Edge**)AC->data)[1-e0->SRMP_WHICH_NODE];
	Energy::Edge* eC = ((Energy::Edge**)AC->data)[e0->SRMP_WHICH_NODE];
	Energy::NonSingletonFactor* A = eA->A;
	Energy::NonSingletonFactor* C = eC->A;
	Energy::Factor* B = eA->B;
	Energy::Edge* e0_rev = (AC->first_out == e0) ? e0->next_out : AC->first_out;
	int a, b, c;
	double delta;
	double* theta = (double*) rbuf.Alloc(B->K*sizeof(double));

	int* TA = (int*) eA->send_message_data;
	int* TAbar = TA + B->K;
	for (b=0; b<B->K; b++)
	{
		double v_min = -e0_rev->m[TA[b]]; // TAbar[a] == 0
		for (a=1; a<A->K/B->K; a++)
		{
			if (v_min > -e0_rev->m[TA[b] + TAbar[a]]) v_min = -e0_rev->m[TA[b] + TAbar[a]];
		}
		theta[b] = v_min;
		if (b==0 || delta>v_min) delta = v_min;
	}
	for (b=0; b<B->K; b++) theta[b] -= delta;

	int* TC = (int*) eC->send_message_data;
	int* TCbar = TC + B->K;
	for (b=0; b<B->K; b++)
	{
		for (c=0; c<C->K/B->K; c++)
		{
			e0->m[TC[b] + TCbar[c]] = theta[b];
		}
	}

	return delta;
}

void PairwiseDualFactorType::SendRestrictedMessage(Energy::Edge* e0)
{
	Energy::NonSingletonFactor* AC = e0->A;
	Energy::Edge* eA = ((Energy::Edge**)AC->data)[1-e0->SRMP_WHICH_NODE];
	Energy::Edge* eC = ((Energy::Edge**)AC->data)[e0->SRMP_WHICH_NODE];
	Energy::NonSingletonFactor* A = eA->A;
	Energy::NonSingletonFactor* C = eC->A;
	Energy::Factor* B = eA->B;
	int a, b, c;

	a = AC->nodes[1-e0->SRMP_WHICH_NODE]->solution;
	b = ComputeRestriction(eA, a);

	int* TC = (int*) eC->send_message_data;
	int* TCbar = TC + B->K;
	for (c=0; c<C->K; c++) e0->m[c] = SRMP_ENERGY_INFTY;
	for (c=0; c<C->K/B->K; c++)
	{
		e0->m[TC[b] + TCbar[c]] = 0;
	}
}


double PairwiseDualFactorType::SendMPLPMessages(Energy::NonSingletonFactor* _A, bool set_solution)
{
	Energy::Edge* eA[2] = { ((Energy::Edge**)_A->data)[0], ((Energy::Edge**)_A->data)[1] };
	Energy::NonSingletonFactor* A[2] = { eA[0]->A, eA[1]->A };
	Energy::Factor* B = eA[0]->B;
	Energy::Node** nodes = _A->nodes;
	int a, b;
	double delta;
	double* theta[2];
	theta[0] = (double*) rbuf.Alloc(2*B->K*sizeof(double));
	theta[1] = theta[0] + B->K;
	int dir;

	for (dir=0; dir<2; dir++)
	{
		int* TA = (int*) eA[dir]->send_message_data;
		int* TAbar = TA + B->K;
		double* rep = nodes[dir]->rep;
		for (b=0; b<B->K; b++)
		{
			double v_min = rep[TA[b]]; // TAbar[a] == 0
			for (a=1; a<A[dir]->K/B->K; a++)
			{
				if (v_min > rep[TA[b] + TAbar[a]]) v_min = rep[TA[b] + TAbar[a]];
			}
			theta[dir][b] = v_min;
		}
	}
	int b_min = 0;
	delta = theta[0][0] + theta[1][0];
	for (b=1; b<B->K; b++)
	{
		if (delta > theta[0][b] + theta[1][b]) { b_min = b; delta = theta[0][b] + theta[1][b]; }
	}
	for (b=0; b<B->K; b++)
	{
		theta[0][b] -= delta;
		theta[1][b] -= delta;
	}

	if (set_solution && (nodes[0]->solution < 0 || nodes[1]->solution < 0))
	{
		if (nodes[0]->solution < 0 && nodes[1]->solution < 0) b = b_min;
		else
		{
			dir = (nodes[0]->solution >= 0) ? 0 : 1;
			b = ComputeRestriction(eA[dir], nodes[dir]->solution);
		}

		for (dir=0; dir<2; dir++)
		if (nodes[dir]->solution < 0)
		{
			int* TA = (int*) eA[dir]->send_message_data;
			int* TAbar = TA + B->K;
			double* rep = nodes[dir]->rep;
			nodes[dir]->solution = TA[b]; // TAbar[0] == 0
			double v_min = rep[TA[b]];
			for (a=1; a<A[dir]->K/B->K; a++)
			{
				if (v_min > rep[TA[b] + TAbar[a]]) { nodes[dir]->solution = TA[b] + TAbar[a]; v_min = rep[TA[b] + TAbar[a]]; }
			}
		}
	}

	int weight[2];
	for (Energy::Edge* e=_A->first_out; e; e=e->next_out) weight[e->SRMP_WHICH_NODE] = e->weight_forward;
	int total_weight = _A->weight_forward + weight[0] + weight[1];
	double total_weight_inv = 1.0 / total_weight;

	for (dir=0; dir<2; dir++)
	{
		double rho = weight[dir]*total_weight_inv;

		int* TA = (int*) eA[dir]->send_message_data;
		int* TAbar = TA + B->K;
		double* rep = nodes[dir]->rep;
		for (b=0; b<B->K; b++)
		{
			for (a=0; a<A[dir]->K/B->K; a++)
			{
				rep[TA[b] + TAbar[a]] = rho*(rep[TA[b] + TAbar[a]] + theta[1-dir][b]);
			}
		}
	}

	return delta;
}









////////////////////////////////////////////////////////////////////////////////

inline int CAST_REP(double* ptr)
{
	return (int)((char*)ptr - (char*)NULL);
}
inline double* CAST_REP(int i)
{
	return (double*)(((char*)NULL) + i);
}

struct PairwiseDualEdgePtr
{
	Energy::Edge* e;
	bool operator<(const PairwiseDualEdgePtr& p) const { return (e->A->rep < p.e->A->rep); }
	bool operator>(const PairwiseDualEdgePtr& p) const { return (e->A->rep > p.e->A->rep); }
};

void Energy::SetFullEdgesDual(int sort_flag)
{
	SetFullEdges(2);

	Node* i;
	NonSingletonFactor* A;
	Factor* B;
	Edge* e;
	dual_sequence = new Sequence;
	Sequence& seq = *dual_sequence;

	seq.num = 0;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (!A->is_removed && !A->first_in) seq.num ++;
	}
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (!i->is_removed && !i->first_in) seq.num ++;
	}

	seq.arr = new FactorPtr[seq.num];
	int t = 0;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (!A->is_removed && !A->first_in) seq.arr[t++].A = A;
	}
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (!i->is_removed && !i->first_in) seq.arr[t++].A = i;
	}

	Options options;
	options.sort_flag = sort_flag;

	SortSequence(seq, options);

	if (!factor_type_general) factor_type_general = new GeneralFactorType;
	for (e=edges->ScanFirst(); e; e=edges->ScanNext())
	{
		e->A->type = factor_type_general;
		e->A->type->InitEdge(e);
	}

	/////////////////////////////////////////////////////////
	factor_type_pairwise_dual = new PairwiseDualFactorType();

	/////////////////////////////////////////////////////////
	// create new Energy instance
	dual_graph = new Energy(seq.num);

	for (t=0; t<seq.num; t++)
	{
		dual_graph->AddNode(seq.arr[t].A->K, seq.arr[t].A->data);
		seq.arr[t].A->rep = CAST_REP(t);
	}

	ReusableBuffer rbuf;

	int phase = -1;
	while ( 1 )
	{
		if (phase<1)
		{
			B = (phase < 0) ? factors.ScanFirst() : factors.ScanNext();
			if (!B) { B = nodes-1; phase = 1; continue; }
			phase = 0;
		}
		else
		{
			B = ((Node*)B) + 1;
			if (B >= nodes+node_num) break;
		}

		if (!B->first_in || !B->first_in->next_in) continue;

		int r, incoming_num;
		for (incoming_num=0, e=B->first_in; e; e=e->next_in) incoming_num ++;
		Energy::Edge** incoming_edges = (Energy::Edge**) rbuf.Alloc(incoming_num*sizeof(Energy::Edge*));
		for (r=0, e=B->first_in; e; e=e->next_in) incoming_edges[r ++] = e;
		quickSort<PairwiseDualEdgePtr>((PairwiseDualEdgePtr*)incoming_edges, 0, incoming_num-1);
		for (r=0; r<incoming_num-1; r++)
		{
			Energy::Edge* e2[2] = { incoming_edges[r], incoming_edges[r+1] };
			Energy::NodeId f2[2] = { CAST_REP(e2[0]->A->rep), CAST_REP(e2[1]->A->rep) };
			dual_graph->AddFactor(2, f2, (double*)e2, factor_type_pairwise_dual);
		}
	}

	dual_graph->primal_graph = this;
}

double Energy::SolveDualGraph(Options& options)
{
	dual_solution_was_inconsistent = false;
	double LB = dual_graph->Solve(options);

	if (!is_cost_best_valid)
	{
		Node* i;
		for (i=dual_graph->nodes; i<dual_graph->nodes+dual_graph->node_num; i++) i->solution = i->solution_best;
		ConvertSolutionDualToPrimal();
	}

	if (options.verbose)
	{
		printf("\n");
		printf("best cost for the orig energy=%f\n", cost_best);
		printf("best cost for the dual energy=%f\n", dual_graph->cost_best);
		if (dual_solution_was_inconsistent)
			printf("  NOTE: there were inconsistent solutions in the dual graph.\n  The cost output above is for the orig energy\n  - see the implementation of Energy::ComputeCost()\n");
	}

	return LB;
}


double Energy::ConvertSolutionDualToPrimal()
{
	int t;
	for (t=0; t<dual_sequence->num; t++)
	{
		Factor* B = dual_sequence->arr[t].A;
		if (B->arity == 1) ((Node*)B)->solution = dual_graph->nodes[t].solution;
		else
		{
			NonSingletonFactor* A = (NonSingletonFactor*) B;
			int r, k = dual_graph->nodes[t].solution;
			for (r=A->arity-1; r>=0; r--)
			{
				int K = A->nodes[r]->K;
				A->nodes[r]->solution = k % K;
				k /= K;
			}
		}
	}

	return ComputeCost();
}

/*
double Energy::ConvertSolutionDualToPrimal()
{
	// copy solution
	Node* i;
	NonSingletonFactor* A;

	for (i=nodes; i<nodes+node_num; i++)
	{
		if (!i->is_removed && !i->first_in) i->solution = dual_graph->nodes[CAST_REP(i->rep)].solution;
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (!A->is_removed && !A->first_in)
		{
			int r, k = dual_graph->nodes[CAST_REP(A->rep)].solution;
			for (r=A->arity-1; r>=0; r--)
			{
				int K = A->nodes[r]->K;
				A->nodes[r]->solution = k % K;
				k /= K;
			}
		}
	}

	return ComputeCost();
}
*/
