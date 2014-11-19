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
#include <srmp/Algs/util.h>
#include <srmp/SRMP.h>
#include <srmp/FactorTypes/GeneralType.h>

using namespace srmpLib;

void Energy::Print()
{
	Node* i;
	NonSingletonFactor* A;
	int k;
	printf("nodes:\n");
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (i->is_removed) continue;
		printf("\t%d :   ", (int)(i-nodes));
		for (k=0; k<i->K; k++) printf("%f ", i->data[k]);
		printf("\n");
	}
	printf("factors:\n");
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		printf("\t");
		for (k=0; k<A->arity; k++)
		{
			printf("%d ", (int)(A->nodes[k]-nodes));
		}
		printf(":  ");
		if (! A->data) printf("-");
		else
		for (k=0; k<A->K; k++) printf("%f ", A->data[k]);
		printf("\n");
	}
}


#ifdef SRMP_TEST_FACTOR_TYPES

void Energy::TestComputePartialReparameterization(NonSingletonFactor* A, double* theta)
{
	A->type->ComputePartialReparameterization(A, theta);

	if (A->type == factor_type_general) return;

	int k;
	double* costs = new double[2*A->K];
	double* theta_test = costs + A->K;
	int* _buf = new int[A->arity];
	GetFactorCosts(A, costs, _buf);

	FactorType* type_old = A->type;
	double* data_old = A->data;

	if (!factor_type_general) factor_type_general = new GeneralFactorType;
	A->type = factor_type_general;
	A->type->InitFactor(A, costs, FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY);
	void** send_message_data_old;
	ReusableBuffer q;
	int n_out = 0;
	Edge* e2;
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		n_out ++;
		send_message_data_old = (void**) q.Realloc(n_out*sizeof(void*));
		send_message_data_old[n_out-1] = e2->send_message_data;
		A->type->InitEdge(e2);
	}

	A->type->ComputePartialReparameterization(A, theta_test);
	for (k=0; k<A->K; k++)
	{
		if (fabs(theta[k] - theta_test[k]) > 1e-6)
		{
			printf("Incorrect output in ComputePartialReparameterization()\n");
			exit(1);
		}
	}

	A->type = type_old;
	A->data = data_old;
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		e2->send_message_data = *send_message_data_old ++;
	}

	delete [] costs;
	delete [] _buf;
}

double Energy::TestSendMessage(Edge* e)
{
	NonSingletonFactor* A = e->A;
	double delta = A->type->SendMessage(e);

	if (A->type == factor_type_general) return delta;

	int k;
	double* costs = new double[A->K + e->B->K];
	int* _buf = new int[A->arity];
	GetFactorCosts(A, costs, _buf);

	FactorType* type_old = A->type;
	double* data_old = A->data;
	double* m_old = e->m;

	e->m = costs + A->K;
	if (!factor_type_general) factor_type_general = new GeneralFactorType;
	A->type = factor_type_general;
	A->type->InitFactor(A, costs, FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY);
	void** send_message_data_old;
	ReusableBuffer q;
	int n_out = 0;
	Edge* e2;
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		n_out ++;
		send_message_data_old = (void**) q.Realloc(n_out*sizeof(void*));
		send_message_data_old[n_out-1] = e2->send_message_data;
		A->type->InitEdge(e2);
	}

	double delta_test = A->type->SendMessage(e);
	if (fabs(delta - delta_test) > 1e-6)
	{
		printf("Incorrect output in SendMessage()\n");
		exit(1);
	}
	for (k=0; k<e->B->K; k++)
	{
		if (fabs(e->m[k] - m_old[k]) > 1e-6)
		{
			printf("Incorrect output in SendMessage()\n");
			exit(1);
		}
	}

	e->m = m_old;
	A->type = type_old;
	A->data = data_old;
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		e2->send_message_data = *send_message_data_old ++;
	}

	delete [] costs;
	delete [] _buf;

	return delta;
}

void Energy::TestSendRestrictedMessage(Edge* e)
{
	NonSingletonFactor* A = e->A;
	A->type->SendRestrictedMessage(e);

	if (A->type == factor_type_general) return;

	int k;
	double* costs = new double[A->K + e->B->K];
	int* _buf = new int[A->arity];
	GetFactorCosts(A, costs, _buf);

	FactorType* type_old = A->type;
	double* data_old = A->data;
	double* m_old = e->m;

	e->m = costs + A->K;
	if (!factor_type_general) factor_type_general = new GeneralFactorType;
	A->type = factor_type_general;
	A->type->InitFactor(A, costs, FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY);
	void** send_message_data_old;
	ReusableBuffer q;
	int n_out = 0;
	Edge* e2;
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		n_out ++;
		send_message_data_old = (void**) q.Realloc(n_out*sizeof(void*));
		send_message_data_old[n_out-1] = e2->send_message_data;
		A->type->InitEdge(e2);
	}

	A->type->SendRestrictedMessage(e);
	int k_opt = 0;
	for (k=1; k<e->B->K; k++) if (e->m[k_opt] > e->m[k]) k_opt = k;
	double delta = e->m[k_opt] - m_old[k_opt];
	for (k=1; k<e->B->K; k++)
	{
		if (fabs(e->m[k] - (m_old[k] + delta)) > 1e-6)
		{
			printf("Incorrect output in SendRestrictedMessage()\n");
			exit(1);
		}
	}

	e->m = m_old;
	A->type = type_old;
	A->data = data_old;
	for (e2=A->first_out; e2; e2=e2->next_out)
	{
		e2->send_message_data = *send_message_data_old ++;
	}

	delete [] costs;
	delete [] _buf;
}

double Energy::TestSendMPLPMessages(NonSingletonFactor* A, bool set_solution)
{
	if (A->type == factor_type_general) return A->type->SendMPLPMessages(A, set_solution);

	Edge* e;

	int buf_size = 2*A->K;
	for (e=A->first_out; e; e=e->next_out) buf_size += e->B->K;

	double* costs = new double[buf_size];
	buf_size = 2*A->K;
	for (e=A->first_out; e; e=e->next_out)
	{
		memcpy(costs+buf_size, e->B->rep, e->B->K*sizeof(double));
		buf_size += e->B->K;
	}

	double delta = A->type->SendMPLPMessages(A, set_solution);

	if (A->rep) memcpy(costs+A->K, A->rep, A->K*sizeof(double));

	buf_size = 2*A->K;
	for (e=A->first_out; e; e=e->next_out)
	{
		// exchange e->B->rep and costs+buf_size
		memcpy(costs, costs+buf_size, e->B->K*sizeof(double));
		memcpy(costs+buf_size, e->B->rep, e->B->K*sizeof(double));
		memcpy(e->B->rep, costs, e->B->K*sizeof(double));
		buf_size += e->B->K;
	}

	/////////////////////////////////////////////

	int* _buf = new int[A->arity];
	GetFactorCosts(A, costs, _buf);

	FactorType* type_old = A->type;
	double* data_old = A->data;

	if (!factor_type_general) factor_type_general = new GeneralFactorType;
	A->type = factor_type_general;
	A->type->InitFactor(A, costs, FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY);
	void** send_message_data_old;
	ReusableBuffer q;
	int n_out = 0;
	for (e=A->first_out; e; e=e->next_out)
	{
		n_out ++;
		send_message_data_old = (void**) q.Realloc(n_out*sizeof(void*));
		send_message_data_old[n_out-1] = e->send_message_data;
		A->type->InitEdge(e);
	}

	double delta_test = A->type->SendMPLPMessages(A, set_solution);
	if (fabs(delta - delta_test) > 1e-6)
	{
		printf("Incorrect output in SendMPLPMessages() [delta]\n");
		exit(1);
	}
	int k;
	buf_size = A->K;
	if (A->rep)
	{
		for (k=0; k<A->K; k++)
		{
			if (fabs(A->rep[k] - costs[k+buf_size]) > 1e-6)
			{
				printf("Incorrect output in SendMPLPMessages() [A->rep]\n");
				exit(1);
			}
		}
	}
	buf_size += A->K;
	for (e=A->first_out; e; e=e->next_out)
	{
		for (k=0; k<e->B->K; k++)
		{
			if (fabs(e->B->rep[k] - costs[k+buf_size]) > 1e-6)
			{
				printf("Incorrect output in SendMPLPMessages() [e->B->rep]\n");
				exit(1);
			}
		}
		buf_size += e->B->K;
	}

	A->type = type_old;
	A->data = data_old;
	for (e=A->first_out; e; e=e->next_out)
	{
		e->send_message_data = *send_message_data_old ++;
	}

	delete [] costs;
	delete [] _buf;

	return delta;
}

#endif



#ifdef SRMP_VNK_DEBUG
void Energy::AddRandomEdges(double prob)
{
	NonSingletonFactor* A;
	Node* i;
	int r, s;
	Block<Factor*>** lists = new Block<Factor*>*[arity_max];
	for (r=0; r<arity_max; r++) lists[r] = new Block<Factor*>(1024);

	for (i=nodes; i<nodes+node_num; i++)
	{
		Factor** ptr = lists[0]->New();
		*ptr = i;
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		Factor** ptr = lists[A->arity-1]->New();
		*ptr = A;
	}
	for (s=1; s<arity_max; s++)
	for (r=0; r<s; r++)
	{
		Factor** Aptr;
		Factor** Bptr;
		for (Aptr=lists[s]->ScanFirst(); Aptr; Aptr=lists[s]->ScanNext())
		{
			for (Bptr=lists[r]->ScanFirst(); Bptr; Bptr=lists[r]->ScanNext())
			{
				if (rand() < prob*RAND_MAX && isSuperset(*Aptr, *Bptr)) AddRelaxationEdge(*Aptr, *Bptr);
			}
		}
	}

	for (r=0; r<arity_max; r++) delete lists[r];

//	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
//	{
//		if (!A->first_out) AddRelaxationEdge(A, A->nodes[0]);
//	}
}

Energy::FactorId Energy::GetFactorId(int arity, NodeId* node_indexes)
{
	if (arity == 1)
	{
		if (node_indexes[0]<0 || node_indexes[0]>=node_num) return NULL;
		return nodes + node_indexes[0];
	}

	int i;
	NonSingletonFactor A;
	NonSingletonFactor* B;

	A.nodes = new Node*[arity];
	A.arity = arity;
	for (i=0; i<arity; i++) A.nodes[i] = nodes + node_indexes[i];
	quickSort<Node*>(A.nodes, 0, arity-1);
	A.is_unsorted = 0;

	for (B=factors.ScanFirst(); B; B=factors.ScanNext())
	{
		if (CompareFactors(&A, B) == 0) return B;
	}
	return NULL;
}
#endif



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/*
#include "FactorTypes/GeneralType.h"
#include "FactorTypes/PairwiseType.h"
#include "FactorTypes/SharedPairwiseType.h"
#include "FactorTypes/PottsType.h"
#include "FactorTypes/PatternType.h"

const int K_min = 2;
const int K_max = 2;

double u_min = 0;
double u_max = 3;
double p_min = 0;
double p_max = 5;

PottsFactorType* potts_debug = new PottsFactorType;
PatternFactorType* pattern_debug = new PatternFactorType;

double Rand(double min, double max)
{
	return (min + ((max-min)*rand())/RAND_MAX);
}
int RandInt(double _min, double _max)
{
	int min = (int)_min, max = (int)_max;
	while ( 1 )
	{
		int x = (int)(min + ((max-min+1.0)*rand())/(RAND_MAX+1.0));
		if (x >= min && x <= max) return x;
	}
}

/////////////////////////////////////////////////////////////////////

void GenerateUnary(int K, double* D)
{
	int k;
	for (k=0; k<K; k++) D[k] = RandInt(u_min, u_max);
}
void GeneratePairwise(int K1, int K2, double* V)
{
	if (K1 == K2)
	{
		double lambda = RandInt(p_min, p_max);
		int ki, kj;
		for (ki=0; ki<K1; ki++)
		for (kj=0; kj<K1; kj++)
		{
			V[K1*ki+kj] = (ki==kj) ? 0 : lambda;
	//		V[K1*ki+kj] = RandInt(p_min, p_max);
		}
	}
	else
	{
		int k;
		for (k=0; k<K1*K2; k++) V[k] = RandInt(p_min, p_max);
	}
}


Energy* GenerateGrid(int sizeX, int sizeY, int high_order=0) // 0: no high-order, 1: zero high-order, 2: non-zero high_order
{
	int x, y, i, j, e, node_num = sizeX*sizeY;
	double buf[K_max*K_max*K_max*K_max];

	Energy* g = new Energy(node_num);
	for (i=0; i<node_num; i++)
	{
		g->AddNode(RandInt(K_min, K_max));
		GenerateUnary(g->GetK(i), buf);
		g->AddUnaryFactor(i, buf);
	}

	for (i=0, y=0; y<sizeY; y++)
	for (x=0; x<sizeX; x++, i++)
	{
		for (e=0; e<2; e++)
		{
			if (e==0) { if (x==sizeX-1) continue; j = i+1; }
			else      { if (y==sizeY-1) continue; j = i+sizeX; }

			GeneratePairwise(g->GetK(i), g->GetK(j), buf);
			g->AddPairwiseFactor(i, j, buf);
		}
		if (high_order && x<sizeX-1 && y<sizeY-1)
		{
			int q[4] = {i,i+1,i+sizeX,i+sizeX+1};
			if (high_order>1)
			{
				int k;
				for (k=0; k<g->GetK(q[0])*g->GetK(q[1])*g->GetK(q[2])*g->GetK(q[3]); k++) buf[k] = RandInt(p_min, p_max);
				g->AddFactor(4, q, buf);
			}
			else g->AddFactor(4, q, NULL);
		}
	}
	return g;
}

Energy* GenerateGrid2(int sizeX, int sizeY) // with triple interactions along scanlines
{
	int x, y, i, j, e, node_num = sizeX*sizeY;
	double buf[K_max*K_max*K_max*K_max];

	Energy* g = new Energy(node_num);
	for (i=0; i<node_num; i++)
	{
		g->AddNode(RandInt(K_min, K_max));
		GenerateUnary(g->GetK(i), buf);
		g->AddUnaryFactor(i, buf);
	}

	for (i=0, y=0; y<sizeY; y++)
	for (x=0; x<sizeX; x++, i++)
	{
		for (e=0; e<2; e++)
		{
			if (e==0) { if (x==sizeX-1) continue; j = i+1; }
			else      { if (y==sizeY-1) continue; j = i+sizeX; }

			GeneratePairwise(g->GetK(i), g->GetK(j), buf);
			g->AddPairwiseFactor(i, j, buf);
		}
		if (x<sizeX-2)
		{
			int q[3] = {i,i+1,i+2};
			if ( 0 )
			{
				int k;
				for (k=0; k<g->GetK(q[0])*g->GetK(q[1])*g->GetK(q[2]); k++) buf[k] = RandInt(p_min, p_max);
				g->AddFactor(3, q, buf);
			}
			else if ( 1 && g->GetK(q[0])==g->GetK(q[1]) && g->GetK(q[0])==g->GetK(q[2]))
			{
				int k = rand() % g->GetK(q[0]);
				int q_labeling[3] = { k, k, k };
				PatternFactorType::Input pattern;
				pattern.cost = -RandInt(0, 3);
				pattern.pattern = q_labeling;
				g->AddFactor(3, q, (double*) &pattern, pattern_debug);
			}
			else g->AddFactor(3, q, NULL);
		}
	}
	return g;
}


Energy* GenerateRandom(int node_num, int factor_num, int arity_max)
{
	ReusableBuffer rbuf;
	int i, j, k;
	int* arr = new int[node_num];
	Energy* g = new Energy(node_num);
	for (i=0; i<node_num; i++) g->AddNode(RandInt(K_min, K_max));
	for ( ; factor_num>0; factor_num--)
	{
		int arity = RandInt(1, arity_max);
		for (i=0; i<arity; i++)
		{
			while ( 1 )
			{
				arr[i] = RandInt(0, node_num-1);
				for (j=0; j<i; j++) if (arr[i] == arr[j]) break;
				if (j == i) break;
			}
		}

		if (arity == 2 && (RandInt(0,10)<6))
		{
			int K0 = g->GetK(arr[0]), K1 = g->GetK(arr[1]);
			if (K0 == K1)
			{
				double lambda = RandInt(0, 10);
				g->AddFactor(arity, arr, &lambda, potts_debug);
				continue;
			}
		}
		int K = 1;
		for (i=0; i<arity; i++) K *= g->GetK(arr[i]);
		double* costs = (double*) rbuf.Alloc(K*sizeof(double));
		for (k=0; k<K; k++) costs[k] = RandInt(0, 10);
		g->AddFactor(arity, arr, (arity == 3) ? 0 : costs);
	}

	delete [] arr;

	return g;
}

int main()
{
	int seed;
	int sizeX = 10, sizeY = 10;


	for (seed=9763; ; seed++)
	{
		printf("seed=%d\n", seed);
		srand(seed);
		//Energy* g = GenerateGrid2(sizeX, sizeY);
		//g->SetMinimalEdges();

		Energy* g = GenerateRandom(5, 10, 5);
		//g->Print();
		//g->AddRandomEdges(Rand(0,1));
		//g->SetFullEdgesDual(0);



		Energy::Options options;
		options.eps = 1e-10;
		options.iter_max = 40;
		//options.sort_flag = 1;

		options.method = Energy::Options::SRMP;
		double LB1 = g->Solve(options);

		options.method = Energy::Options::MPLP;
		double LB2 = g->Solve(options);

		options.method = Energy::Options::CMP;
		double LB3 = g->Solve(options);

		delete g;
	}

	return 0;
}
*/
