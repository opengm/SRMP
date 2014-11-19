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
#include "PQ.h"

using namespace srmpLib;

int Energy::ScoreNodeOrdering()
{
	NonSingletonFactor* A;
	int score = 0;
	int i;

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		int min = A->nodes[0]->solution, max = A->nodes[0]->solution;
		for (i=1; i<A->arity; i++)
		{
			if (min > A->nodes[i]->solution) min = A->nodes[i]->solution;
			if (max < A->nodes[i]->solution) max = A->nodes[i]->solution;
		}
		score += max - min;
	}

	return score;
}

struct FactorListItem
{
	Energy::NonSingletonFactor* A;
	FactorListItem* next;
};

Energy::Node* Energy::_GenerateNodeOrdering(Node* seed)
{
	int i, k, assigned_num = 0, counter = 0;
	PriorityQueue<double> pq;
	void* pq_buf = pq.AllocateBuf();
	pq.Reset();
	PriorityQueue<double>::Item* items = new PriorityQueue<double>::Item[node_num];
	PriorityQueue<double>::Item* t;
	FactorListItem* f;

	// Item::slack:
	// < 0 - in the queue
	// > 1 - removed from the queue
	// = 0 - never been in the queue

	for (i=0; i<node_num; i++)
	{
		items[i].slack = 0;
		nodes[i].solution = -1;
	}
	i = (int)(seed - nodes);
	items[i].slack = -1;
	pq.Add(&items[i]);
	while ( 1 )
	{
		t = pq.GetMin();
		if (!t)
		{
			while (items[counter].slack != 0) counter ++;
			items[counter].slack = -1;
			pq.Add(&items[counter]);
			continue;
		}
		pq.SRMP_Remove(t, pq_buf);
		t->slack = 1;
		Node* p = nodes + (int)(t - items);
		//printf("%d ", (int)(t-items));
		p->solution = assigned_num ++;
		if (assigned_num == node_num)
		{
			return p;
		}
		for (f=(FactorListItem*)p->rep; f; f=f->next)
		{
			for (k=0; k<f->A->arity; k++)
			{
				if (f->A->nodes[k]->solution >= 0) continue;
				i = (int)(f->A->nodes[k] - nodes);
				if (items[i].slack == 0)
				{
					items[i].slack = -2 + (1.0 + assigned_num) / node_num ; 
					pq.Add(&items[i]);
				}
				else
				{
					items[i].slack --;
					pq.Decrease(&items[i], &items[i], pq_buf);
				}
			}
		}
	}
}

void generate_permutation(int *buf, int n)
{
	int i, j;

	for (i=0; i<n; i++) buf[i] = i;
	for (i=0; i<n-1; i++)
	{
		do
		{
			j = i + (int) (((double)rand()/(RAND_MAX+1.0))*(n - i));
		} while (j<i || j>=n);
		int tmp = buf[i]; buf[i] = buf[j]; buf[j] = tmp;
	}
}

void Energy::GenerateNodeOrdering(Options& options)
{
	if (node_num == 0) return;

	Node* i;
	NonSingletonFactor* A;
	int k;

	int* order_best = new int[node_num];
	for (k=0; k<node_num; k++) order_best[k] = nodes[k].solution = k;
	int score_best = ScoreNodeOrdering();
	if (options.verbose) printf("orig ordering score: %d\n", score_best);

	if (options.sort_flag > 1)
	{
		srand((unsigned) options.sort_flag);
		generate_permutation(order_best, node_num);
		int score = ScoreNodeOrdering();
		if (options.verbose) printf("random ordering: %d\n", score);
		for (k=0; k<node_num; k++) nodes[k].solution = order_best[k];
		delete [] order_best;
		return;
	}

	Block<FactorListItem> factor_list_items(1024);

	for (i=nodes; i<nodes+node_num; i++) i->rep = NULL;

	for (A=factors.ScanFirst(); A; A=factors.ScanNext())
	{
		for (k=0; k<A->arity; k++)
		{
			i = A->nodes[k];
			FactorListItem* f = factor_list_items.New();
			f->A = A;
			f->next = (FactorListItem*) i->rep;
			i->rep = (double*) f;
		}
	}

	int iter;
	const int ITER_NUM = 5;
	Node* seeds[ITER_NUM];
	seeds[0] = &nodes[0];
	for (iter=0; ; iter++)
	{
		if (options.verbose) printf("greedy from node %d: ", (int)(seeds[iter] - nodes));
		i = _GenerateNodeOrdering(seeds[iter]);
		int score = ScoreNodeOrdering();
		if (options.verbose) printf("%d\n", score);
		if (score_best > score)
		{
			score_best = score;
			for (k=0; k<node_num; k++) order_best[k] = nodes[k].solution;
		}
		if (iter == ITER_NUM - 1) break;
		for (k=0; k<=iter; k++)
		{
			if (i == seeds[k]) break;
		}
		if (k <= iter) break;
		seeds[iter+1] = i;
	}

	for (k=0; k<node_num; k++) nodes[k].solution = order_best[k];
	delete [] order_best;
}

void Energy::SortSequence(Sequence& seq, Options& options)
{
	if (options.sort_flag < 0) return;
	if (options.sort_flag == 0)
	{
		quickSort<FactorPtr>(seq.arr, 0, seq.num-1);
		return;
	}
	GenerateNodeOrdering(options);
	quickSort<FactorPtrX>((FactorPtrX*)seq.arr, 0, seq.num-1);
}

