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

using namespace srmpLib;

/*
struct MyNode;
struct MyEdge
{
	int i; // if edge = { i0, i1 } then i=max{i0,i1}
	Energy::NonSingletonFactor* edge;

	bool operator<(const MyEdge& B) const { return (i < B.i); }
	bool operator>(const MyEdge& B) const { return (i > B.i); }
};

struct MyNode
{
	int neighbor_num; // only higher-order neihbors are counted
	MyEdge* neighbors; // of size neighbor_num, in the decreasing order
	Energy::NonSingletonFactor* flag;
};



void Energy::FindTriplets()
{
	if (triplets.ScanFirst()) return;

	MyNode* p;
	MyNode* q;
	Node* i;
	Node* j;
	Edge* e;
	int k, k2, n;

	Buffer my_buf(1024);
	ReusableBuffer my_rbuf;
	MyNode* my_nodes = new MyNode[node_num];

	for (p=my_nodes, i=nodes; i<nodes+node_num; p++, i++)
	{
		MyEdge* neighbors;
		n = 0;
		for (e=i->first_in; e; e=e->next_in)
		{
			if (e->A->arity != 2) continue;
			if (!e->A->type->incoming_edges_allowed) continue;
			Edge* e2 = (e->A->first_out == e) ? e->next_out : e->A->first_out;
			if (!e2) continue;
			j = (Node*)e2->B;
			if (j < i) continue;
			neighbors = (MyEdge*) my_rbuf.Realloc((n+1)*sizeof(MyEdge));
			neighbors[n].i = (int)(j - nodes);
			neighbors[n].edge = e->A;
			n ++;
		}
		if (n > 0) quickSort<MyEdge>(neighbors, 0, n-1);
		p->neighbor_num = n;
		p->neighbors = (MyEdge*) my_buf.Alloc(n*sizeof(MyEdge));
		for (k=0; k<n; k++) p->neighbors[k] = neighbors[n-1-k];

		p->flag = NULL;
	}

	for (p=my_nodes; p<my_nodes+node_num; p++)
	{
		for (k=0; k<p->neighbor_num; k++)
		{
			q = &my_nodes[p->neighbors[k].i];
			q->flag = p->neighbors[k].edge;
			for (k2=0; k2<q->neighbor_num; k2++)
			{
				if (my_nodes[q->neighbors[k2].i].flag)
				{
					// found a triplet
					Triplet* t = triplets.New();
					t->e[0] = p->neighbors[k].edge;
					t->e[1] = q->neighbors[k2].edge;
					t->e[2] = my_nodes[q->neighbors[k2].i].flag;
					t->gap = 0;
				}
			}
		}
		for (k=0; k<p->neighbor_num; k++) my_nodes[p->neighbors[k].i].flag = NULL;
	}

	delete [] my_nodes;
}

void Energy::Tighten()
{
	FindTriplets();

	Node* i;
	Node* j;
	Node* k;
	NonSingletonFactor* A;
	Edge* e;
	Buffer my_buf(512);
	ReusableBuffer my_rbuf;
	double gap;
	int a;
	Triplet* t;

	// compute LPDG, write it to Factor::rep
	for (i=nodes; i<nodes+node_num; i++)
	{
		i->solution = (is_solution_best_initialized) ? i->solution_best : 0;
		gap = i->data[0];
		for (a=1; a<i->K; a++) { if (gap > i->data[a]) gap = i->data[a]; }
		gap -= i->data[i->solution];
		i->rep = (double*)my_buf.Alloc(sizeof(double));
		*i->rep = gap;
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->arity != 2) continue;
		if (!A->type->incoming_edges_allowed) continue;
		double* theta = (double*) my_rbuf.Alloc(A->K*sizeof(double));
		SRMP_COMPUTE_PARTIAL_REPARAMETERIZATION(A, theta);
		for (e=A->first_in; e; e=e->next_in)
		{
			for (a=0; a<A->K; a++) theta[a] += e->m[a];
		}
		gap = theta[0];
		for (a=1; a<A->K; a++) { if (gap > theta[a]) gap = theta[a]; }
		gap -= theta[A->nodes[0]->solution*A->nodes[0]->K + A->nodes[1]->solution];
		i->rep = (double*)my_buf.Alloc(sizeof(double));
		*i->rep = gap;
	}

	/////////////////////////////////////////////////////////
	int total_num = 0, nonzero_num = 0;
	for (t=triplets.ScanFirst(); t; t=triplets.ScanNext())
	{
		if (t->gap < 0) continue; // this triplet has already been added
		total_num ++;
		t->GetNodes(i, j, k);
		t->gap = *i->rep + *j->rep + *k->rep + *t->e[0]->rep + *t->e[1]->rep + *t->e[2]->rep;
		if (t->gap > 0) nonzero_num ++;
	}
	printf("%d triplets (%d non-zero)\n", total_num, nonzero_num);
}
*/
