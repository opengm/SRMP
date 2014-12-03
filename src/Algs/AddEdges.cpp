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
#include <srmp/FactorTypes/PairwiseType.h>
#include <srmp/FactorTypes/GeneralType.h>


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

using namespace srmpLib;

inline void Energy::IntersectFactors(NonSingletonFactor* A, NonSingletonFactor* B, NonSingletonFactor* C)
{
	C->arity = 0;
	C->is_unsorted = 0;
	C->is_removed = 0;
	int i=0, j=0;
	Node** Anodes = A->nodes + A->is_unsorted*A->arity;
	Node** Bnodes = B->nodes + B->is_unsorted*B->arity;
	while ( i<A->arity && j<B->arity )
	{
		if (Anodes[i] == Bnodes[j]) { C->nodes[C->arity ++] = Anodes[i]; i++; j++; }
		else if (Anodes[i] < Bnodes[j]) i++;
		else j++;
	}
}

inline bool Energy::isSuperset(NonSingletonFactor* A, NonSingletonFactor* B)
{
	int i=0, j=0;
	Node** Anodes = A->nodes + A->is_unsorted*A->arity;
	Node** Bnodes = B->nodes + B->is_unsorted*B->arity;
	while ( i<A->arity && j<B->arity )
	{
		if (Anodes[i]  > Bnodes[j]) return false;
		if (Anodes[i] == Bnodes[j]) j++;
		i++;
	}
	return (j==B->arity) ? true : false;
}

bool Energy::isSuperset(FactorId _A, FactorId _B)
{
	NonSingletonFactor* A = (NonSingletonFactor*)_A;
	Factor* B = (Factor*)_B;

	if (A->arity <= B->arity) return false;

	int i, j;
	Node** Anodes_sorted = GetSortedNodesPtr(A);
	Node** Bnodes_sorted = GetSortedNodesPtr(B);
	for (i=j=0; i<A->arity && j<B->arity; i++)
	{
		if (Anodes_sorted[i] > Bnodes_sorted[j]) return false;
		if (Anodes_sorted[i] == Bnodes_sorted[j]) j ++;
	}
	if (j < B->arity) return false;

	return true;
}

void Energy::AddRelaxationEdge(FactorId _A, FactorId _B)
{
	AllocateEdges();

	NonSingletonFactor* A = (NonSingletonFactor*)_A;
	Factor* B = (Factor*)_B;
	Edge* e;

	if (A->arity <= 3)
	{
		// check whether the edge is already present
		for (e=A->first_out; e; e=e->next_out)
		{
			if (e->B == B) return;
		}
	}

	e = edges->New();
	e->A = A;
	e->B = B;
	e->next_out = A->first_out;
	A->first_out = e;
	e->next_in = B->first_in;
	B->first_in = e;
	e->m = NULL;
}

void Energy::SetMinimalEdges()
{
	if (edges)
	{
		printf("SetMinimalEdges() cannot be called: edges already exist\n");
		exit(1);
	}
	AllocateEdges();

	int k;
	NonSingletonFactor* A;
	Node* i;

	for (i=nodes; i<nodes+node_num; i++) i->first_in = NULL;

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		assert(A->arity > 1);

		A->first_out = NULL;
		for (k=0; k<A->arity; k++)
		{
			i = A->nodes[k];

			Edge* e = edges->New();
			e->A = A;
			e->B = i;
			e->next_out = A->first_out;
			A->first_out = e;
			e->next_in = i->first_in;
			i->first_in = e;
			e->m = NULL;
		}
	}
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

struct FactorListItem
{
	Energy::NonSingletonFactor* A;
	FactorListItem* next;
};

class NonSingletonFactorOrderedList
{
public:
	NonSingletonFactorOrderedList(int arity_max);
	~NonSingletonFactorOrderedList();

	void Add(Energy::NonSingletonFactor* A); // use A->rep for maintaining the list
	Energy::NonSingletonFactor* RemoveLargest();

private:
	int arity_max, r_current;
	Energy::NonSingletonFactor** lists_begin;
	Energy::NonSingletonFactor** lists_end;
};

NonSingletonFactorOrderedList::NonSingletonFactorOrderedList(int _arity_max)
	: arity_max(_arity_max), r_current(_arity_max)
{
	int r;
	lists_begin = new Energy::NonSingletonFactor*[2*(arity_max+1)];
	lists_end = lists_begin + (arity_max+1);
	for (r=1; r<=arity_max; r++) lists_begin[r] = lists_end[r] = NULL;
}
NonSingletonFactorOrderedList::~NonSingletonFactorOrderedList()
{
	delete [] lists_begin;
}
void NonSingletonFactorOrderedList::Add(Energy::NonSingletonFactor* A)
{
	if (lists_end[A->arity]) lists_end[A->arity]->rep = (double*)A;
	else lists_begin[A->arity] = A;
	lists_end[A->arity] = A;
}
Energy::NonSingletonFactor* NonSingletonFactorOrderedList::RemoveLargest()
{
	while ( 1 )
	{
		if (lists_begin[r_current])
		{
			Energy::NonSingletonFactor* A = lists_begin[r_current];
			lists_begin[r_current] = (lists_begin[r_current] == lists_end[r_current]) ? 
				NULL : (Energy::NonSingletonFactor*) A->rep;
			return A;
		}
		r_current --;
		if (r_current == 0) return NULL;
	}
}

// idea:
// For each singleton factor {i} keep the list of factors that contain in.
// Initialize with all such factors, but then remove factors A for which
// there exists an intermediate factor B with {v} \subset B \subset A.
// To this, go through non-singleton factors A starting from the largest and do the following: 
// - compute non-singleton intersections of A with other factors (using lists stored at singleton factors)
// - add them with zero potentials, if they don't exist yet
// - select those intersections which are maximal 
// - for each such maximal intersection B \subset A add edge A->B and remove A from the list of each singleton factor {i}\subset B
//
// As the final step, go through singleton factors and add edges A->{i} for all A in the list for {i}

void Energy::SetMaximalEdges()
{
	int t, r, rA;
	NonSingletonFactor* A;
	NonSingletonFactor* B;
	NonSingletonFactor C; // intersection of A and B
	Node* i;
	Edge* e;
	NonSingletonFactorOrderedList L(arity_max);
	FactorListItem* p;
	FactorListItem* q;
	DBlock<FactorListItem> factor_list_items(512);
	FactorListItem** children;

	C.nodes = new Node*[arity_max];
	children = new FactorListItem*[arity_max];

	// use NonSingletonFactor::first_out for a list of factors (of a given arity); cast types between (Edge*) and (NonSingletonFactor*)
	// use Node::first_in for keeping a list of FactorListItem's that contain this node; cast types between (Edge*) and (FactorListItem*)

	// set Node::tmp_list's
	for (i=nodes; i<nodes+node_num; i++)
	{
		i->tmp1 = 0;
		i->first_in = NULL;
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		A->first_in = A->first_out = NULL;
		if (A->is_removed) continue;

		if (A->arity > 2) L.Add(A);

		for (r=0; r<A->arity; r++)
		{
			p = factor_list_items.New();
			p->A = A;
			p->next = (FactorListItem*)A->nodes[r]->first_in;
			A->nodes[r]->first_in = (Edge*)p;
		}
//PrintFactor(A); printf(" ");
	}
//printf(" | ");
	// now go through non-singleton factors starting from the largest
	while ((A=L.RemoveLargest()))
	{
		rA = A->arity;

		for (r=2; r<rA; r++) children[r] = NULL;
		// first, set children to be the list of proper non-singleton subsets of A that already exist
		for (r=0; r<rA; r++)
		for (p=(FactorListItem*)A->nodes[r]->first_in; p; p=p->next)
		{
			B = p->A;
			if (A->arity == B->arity || B->arity==1) continue;
			if (!isSuperset(A, B)) continue;
			if (A->nodes[r] != B->nodes[B->is_unsorted*B->arity]) continue;

			q = factor_list_items.New();
			q->A = B;
			q->next = children[B->arity];
			children[B->arity] = q;
		}
		// now add new children (overlaps with other factors)
		for (r=0; r<rA; r++)
		for (p=(FactorListItem*)A->nodes[r]->first_in; p; p=p->next)
		{
			B = p->A;
			if (A == B) continue;
			IntersectFactors(A, B, &C);
			if (C.arity == 1 || C.arity == A->arity) continue;

			// check whether C is already present
			for (q=children[C.arity]; q; q=q->next) if (CompareFactors(q->A, &C) == 0) break;
			if (q) continue;

			NonSingletonFactor* _C =  factors.New();
			_C->arity = C.arity;
			_C->is_unsorted = 0;
			_C->is_removed = 0;
			_C->K = 1;
			_C->nodes = (Node**) buf.Alloc(C.arity*sizeof(Node*));
			for (t=0; t<C.arity; t++)
			{
				_C->nodes[t] = C.nodes[t];
				_C->K *= C.nodes[t]->K;

				q = factor_list_items.New();
				q->A = _C;
				q->next = (FactorListItem*)C.nodes[t]->first_in;
				C.nodes[t]->first_in = (Edge*)q;
			}
			_C->first_in = _C->first_out = NULL;
			InitFactor(_C, NULL);

			if (_C->arity > 2)
			{
				L.Add(_C);
			}

//PrintFactor(_C); printf(" ");

			q = factor_list_items.New();
			q->A = _C;
			q->next = children[_C->arity];
			children[_C->arity] = q;
		}
		// select maximal factors in 'children'
		A->first_out = NULL;
		for (r=rA-1; r>1; r--)
		for (p=children[r]; p; p=q)
		{
			for (e=A->first_out; e; e=e->next_out)
			{
				if (e->B->arity>1 && isSuperset((NonSingletonFactor*)e->B, p->A)) break;
			}
			if (!e)
			{
				e = edges->New();
				e->A = A;
				e->B = p->A;
				e->next_out = e->A->first_out;
				e->A->first_out = e;
				e->next_in = e->B->first_in;
				e->B->first_in = e;
				e->m = NULL;

//PrintEdge(e); printf(" ");

				for (t=0; t<p->A->arity; t++) p->A->nodes[t]->tmp1 = 1; // mark such nodes
			}
			q = p->next;
			factor_list_items.Delete(p);
		}
		// remove A from the list of appropriate singleton factors
		for (r=0; r<rA; r++)
		{
			i = A->nodes[r];
			if (!i->tmp1) continue;
			i->tmp1 = 0;
			FactorListItem** ptr = (FactorListItem**)&i->first_in;
			while ((*ptr)->A != A) ptr = &(*ptr)->next;
			*ptr = (*ptr)->next;
		}
	}

	for (i=nodes; i<nodes+node_num; i++)
	{
		for (p=(FactorListItem*)i->first_in, i->first_in=NULL; p; p=p->next)
		{
			e = edges->New();
			e->A = p->A;
			e->B = i;
			e->next_out = e->A->first_out;
			e->A->first_out = e;
			e->next_in = e->B->first_in;
			e->B->first_in = e;
			e->m = NULL;

//PrintEdge(e); printf(" ");
		}
	}
//printf("\n");

	delete [] C.nodes;
	delete [] children;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// A and B must be the same factor (possibly, with a different order of nodes).
// A and B must be of a general type with non-zero costs.
// t must be of 'GeneralFactorType'
// The function moves costs of B to A
void MergeFactorCosts(Energy::NonSingletonFactor* A, Energy::NonSingletonFactor* B, GeneralFactorType* t)
{
	int k;
	Energy::NonSingletonFactor _A;
	Energy::Edge e;

	_A.arity = A->arity;
	_A.K = A->K;
	_A.nodes = A->nodes;
	_A.is_unsorted = A->is_unsorted;
	_A.data = A->data;
	_A.type = t;

	_A.first_in = NULL;
	_A.first_out = &e;

	e.next_out = NULL;
	e.A = A;
	e.B = B;
	t->InitEdge(&e);

	e.m = B->data;
	for (k=0; k<A->K; k++) e.m[k] = -e.m[k];
	t->ComputePartialReparameterization(&_A, _A.data);
}

void Energy::MoveCostsToOuterFactors()
{
	if (edges)
	{
		printf("MoveCostsToOuterFactors() cannot be called: edges already exist\n");
		exit(1);
	}
	AllocateEdges();

	int r, r2, rA, k;
	Node* i;
	Edge* e;
	FactorListItem* p;
	NonSingletonFactor* A;
	NonSingletonFactorOrderedList L(arity_max);
	DBlock<FactorListItem> factor_list_items(512);
	ReusableBuffer rbuf;
	Buffer buf_local(1024);
	GeneralFactorType general_tmp;

	for (i=nodes; i<nodes+node_num; i++) i->tmp1 = 0;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		// convert to general
		if (A->type == factor_type_pairwise || A->type == factor_type_general)
		{
			if (!A->data)
			{
				A->data = (double*) buf.Alloc(A->K*sizeof(double));
				memset(A->data, 0, A->K*sizeof(double));
			}
		}
		else
		{
			double* costs = (double*) rbuf.Alloc(A->K*sizeof(double) + A->arity*sizeof(int));
			GetFactorCosts(A, costs, (void*)(costs + A->K));
			InitFactor(A, costs, NULL);
		}

		A->first_in = A->first_out = NULL;

		L.Add(A);
	}

	while ((A=L.RemoveLargest()))
	{
		rA = A->arity;

		// find factors containing A
		for (r=0; r<rA; r++) A->nodes[r]->tmp1 = 1;
		for (r=0; r<rA; r++)
		{
			for (p=(FactorListItem*)A->nodes[r]->first_in; p; p=p->next)
			{
				if (p->A == A) continue;
				for (r2=0; !p->A->nodes[r2]->tmp1; r2++) {}
				if (p->A->nodes[r2] != A->nodes[r]) continue;

				if (isSuperset(p->A, A))
				{
					if (p->A->arity > A->arity) AddRelaxationEdge(p->A, A);
					else MergeFactorCosts(p->A, A, &general_tmp);
					A->is_removed = 1;
				}
			}
		}
		for (r=0; r<rA; r++) A->nodes[r]->tmp1 = 0;

		if (A->is_removed)
		{
			if (A->first_in)
			{
				int num_in = 0;
				for (e=A->first_in; e; e=e->next_in) num_in ++;
				for (k=0; k<A->K; k++) A->data[k] /= -num_in;
				for (e=A->first_in; e; e=e->next_in)
				{
					e->A->type->InitEdge(e);
					e->m = (double*) buf_local.Alloc(A->K*sizeof(double));
					memcpy(e->m, A->data, A->K*sizeof(double));
				}
			}
			memset(A->data, 0, A->K*sizeof(double));

			continue;
		}

		for (r=0; r<rA; r++)
		{
			p = factor_list_items.New();
			p->A = A;
			p->next = (FactorListItem*)A->nodes[r]->first_in;
			A->nodes[r]->first_in = (Edge*)p;
		}
	}
	for (i=nodes; i<nodes+node_num; i++)
	{
		p = (FactorListItem*)i->first_in;
		i->first_in = NULL;
		if (!p) continue;

		int num_in = 0;
		for ( ; p; p=p->next)
		{
			AddRelaxationEdge(p->A, i);
			num_in ++;
		}
		for (k=0; k<i->K; k++) i->data[k] /= -num_in;
		for (e=i->first_in; e; e=e->next_in)
		{
			e->A->type->InitEdge(e);
			e->m = (double*) buf_local.Alloc(i->K*sizeof(double));
			memcpy(e->m, i->data, i->K*sizeof(double));
		}
		memset(i->data, 0, i->K*sizeof(double));
		i->is_removed = (num_in > 1) ? 1 : 0; // if i belongs to a single outer factor then don't remove it, so that the solution for it is set correctly
	}

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		A->type->ComputePartialReparameterization(A, A->data);
	}

	delete edges;
	edges = NULL;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

class UnionFind // uses Factor::rep
{
public:
	void Add(Energy::Factor* A);
	void Add(Energy::Factor* A, Energy::Factor* B); // add A to the same component as B
	bool CheckAndMerge(Energy::Factor* A, Energy::Factor* B); // returns true if A, B were already in the same component
};

inline void UnionFind::Add(Energy::Factor* A)
{
	A->rep = (double*) A;
}

inline void UnionFind::Add(Energy::Factor* A, Energy::Factor* B)
{
	A->rep = (double*) B;
}

bool UnionFind::CheckAndMerge(Energy::Factor* A, Energy::Factor* B)
{
	Energy::Factor* A_parent = A;
	while (A_parent != (Energy::Factor*)A_parent->rep) A_parent = (Energy::Factor*)A_parent->rep;
	Energy::Factor* B_parent = B;
	while (B_parent != (Energy::Factor*)B_parent->rep) B_parent = (Energy::Factor*)B_parent->rep;

	bool result;
	if (A_parent == B_parent) result = true;
	else
	{
		result = false;
		B_parent->rep = (double*) A_parent;
	}

	// path compression
	while (A != A_parent)
	{
		Energy::Factor* A_next = (Energy::Factor*) A->rep;
		A->rep = (double*) A_parent;
		A = A_next;
	}
	A = B;
	while (A != A_parent)
	{
		Energy::Factor* A_next = (Energy::Factor*) A->rep;
		A->rep = (double*) A_parent;
		A = A_next;
	}

	return result;
}

struct NestedEdgePair
{
	Energy::NonSingletonFactor* A;
	Energy::Factor* B;
};

void Energy::SetFullEdges(int method)
{
	if (edges)
	{
		printf("SetFullEdges() cannot be called: edges already exist\n");
		exit(1);
	}

	if (method > 0) MoveCostsToOuterFactors();

	AllocateEdges();
	SetMaximalEdges();
	if (method==0) return;

	NonSingletonFactorOrderedList L(arity_max);
	UnionFind U;
	NonSingletonFactor* A;
	Factor* B = NULL;
	Edge* e;
	int r;

	//for (e=edges->ScanFirst(); e; e=edges->ScanNext()) { PrintEdge(e); printf("\n"); }

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed) continue;
		L.Add(A);
	}

	// remove unnecessary edges

	Block<NestedEdgePair> new_edges(512);

	r = arity_max + 1;
	while ( 1 )
	{
		if (r > 1)
		{
			B = L.RemoveLargest();
			if (B) r = B->arity;
			else { r = 1; B = nodes-1; continue; }
		}
		else
		{
			B = ((Node*)B) + 1;
			if (B >= nodes+node_num) break;
		}

		if (!B->first_in)  // outer factor
		{
			B->rep = NULL;
			continue;
		}

		// traverse all parents of B, add them to U
		A = B->first_in->A;
		U.Add(A);
		NonSingletonFactor* list_begin = A;
		NonSingletonFactor* list_end = A;
		for (e=B->first_in->next_in; e; e=e->next_in)
		{
			A = e->A;
			U.Add(A);
			list_end->first_out = (Edge*) A;
			list_end = A;
		}
		while ((A=list_begin))
		{
			list_begin = (A == list_end) ? NULL : (NonSingletonFactor*)A->first_out;
			for (e=A->first_in; e; e=e->next_in)
			{
				if (e->A->rep) // already in the list
				{
					U.CheckAndMerge(A, e->A);
				}
				else
				{
					U.Add(e->A, A);
					list_end->first_out = (Edge*) e->A; 
					list_end = e->A;
					if (!list_begin) list_begin = e->A;
				}
			}
		}

		// now analyze connected components
		int component_num = 1;
		e = B->first_in;
		U.Add(B, e->A);
		for (e=e->next_in; e; e=e->next_in)
		{
			if (U.CheckAndMerge(B, e->A)) continue;

			component_num ++;
			NestedEdgePair* q = new_edges.New();
			q->A = e->A; q->B = B;
			if (method == 2)
			{
				while (q->A->first_in) q->A = q->A->first_in->A;
			}
		}
		if (component_num > 1)
		{
			e = B->first_in;

			NestedEdgePair* q = new_edges.New();
			q->A = e->A; q->B = B;
			if (method == 2)
			{
				while (q->A->first_in) q->A = q->A->first_in->A;
			}
		}

		// finally, set 'rep' pointers back to NULL
		B->rep = NULL;
		A = B->first_in->A;
		while ( 1 )
		{
			A->rep = NULL;
			if (A == list_end) break;
			A = (NonSingletonFactor*) A->first_out;
		}
	}

	delete edges;
	edges = NULL;
	AllocateEdges();

	Node* i;
	for (i=nodes; i<nodes+node_num; i++) i->first_in = NULL;
	for (A=factors.ScanFirst(); A; A=factors.ScanNext()) A->first_in = A->first_out = NULL;

	NestedEdgePair* q;
	for (q=new_edges.ScanFirst(); q; q=new_edges.ScanNext())
	{
		AddRelaxationEdge(q->A, q->B);
	}
	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		if (A->is_removed || A->first_in) continue;
		{
			for (r=0; r<A->arity; r++)
			{
				if (!A->nodes[r]->is_removed) AddRelaxationEdge(A, A->nodes[r]); // this node is contained only in A. Add this edge so that the solution for this node is assigned correctly
			}
		}
	}

	//printf("new edges:\n");
	//for (e=edges->ScanFirst(); e; e=edges->ScanNext()) { PrintEdge(e); printf("\n"); }
}

