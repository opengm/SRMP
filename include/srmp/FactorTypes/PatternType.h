/* PatternType.h */
/* Copyright Vladimir Kolmogorov vnk@ist.ac.at */
/* Created June 2013 */

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

/********************************************************************
 Example of a factor with incremental updates. In other words, it exploits the following fact:
 when the code calls
 PrepareFactor(A), SendMessage(e_1), SendMessage(e_2), ..., SendMessage(e_r), PrepareFactor(A), ...
 then between calls SendMessage(e_k) and SendMessage(e_{k+1}) only message e_k->m changes.

 Note: Ideally, we also would like to do incremental updates for SendRestrictedMessage() [which is used for extracting a primal solution].
 However, the current code structure doesn't support that. (Maybe, some time in the future...)
 Therefore, primal solution computations should not be called every 3 iterations.
********************************************************************/







/* Example usage: minimize f(x,y,z) = 2*x + x*y + P(x,y,z)
   where x,y,z \in {0,1,2} and P(x,y,z)=-5 if (x,y,z)=(1,0,1), and P(x,y,z)=0 otherwise.

#include <stdio.h>
#include "SRMP.h"
#include "PatternType.h"

int main(int argc, char* argv[])
{
	int i, a, b, node_num = 3;
	int f[3];
	int pattern[3];
	double lambda, D[3], V[3*3];

	Energy* g = new Energy(node_num);
	PatternFactorType* pattern_type = new PatternFactorType;

	for (i=0; i<node_num; i++)
	{
		g->AddNode(3); // add node with 3 labels
	}

	// add unary term 2*x
	D[0] = 0; D[1] = 2; D[2] = 4;
	g->AddUnaryFactor(0, D);

	// add pairwise term x*y
	for (a=0; a<3; a++)
	for (b=0; b<3; b++)
	{
		V[3*a + b] = a*b;
	}
	g->AddPairwiseFactor(0, 1, V);

	// add term P(x,y,z)
	PatternFactorType::Input input;
	f[0] = 0; f[1] = 1; f[2] = 2;
	pattern[0] = 1; pattern[1] = 0; pattern[2] = 1;
	input.cost = -5; // must be non-positive
	input.pattern = pattern;
	g->AddFactor(2, f, (double*)&input, pattern_type);

	// call solver
	Energy::Options options;
	options.iter_max = 10;
	g->Solve(options);

	// print solution
	for (i=0; i<node_num; i++)
	{
		printf("%d ", g->GetSolution(i));
	}

	// done
	delete g;
	delete pattern_type;

	return 0;
}

*/


#ifndef SRMP_FUAHNGIUAGKAFASJGASG
#define SRMP_FUAHNGIUAGKAFASJGASG

#include "../SRMP.h"

namespace srmpLib {

struct PatternFactorType : Energy::FactorType
{
	struct Input
	{
		double cost;
		int* pattern; // array of size 'A->arity'
	};

	PatternFactorType();
	~PatternFactorType();

	// user_data points to 'Input' structure
	void InitFactor(Energy::NonSingletonFactor* A, double* user_data, unsigned flags);
	double GetCost(Energy::NonSingletonFactor* A);

	void ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* rep);

	void InitEdge(Energy::Edge* e);

	double SendMessage(Energy::Edge* e);
	void SendRestrictedMessage(Energy::Edge* e);

	double SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution);

	bool PrepareFactor(Energy::NonSingletonFactor* A);

	/////////////////////////////////////////////////////
private:
	Buffer buf;
	ReusableBuffer rbuf;

	// A->data actually points to a 'FactorData' structure
	struct FactorData
	{
		double cost; // <= 0
		int* pattern; // array of labels of size 'arity'

		double sum_min;     // =        sum_{i != index_last} \min_k { -m_i[k] }
		double sum_pattern; // = cost + sum_{i != index_last} { -m_i[pattern[i]] }
		Energy::Edge* e_last; // index_last = (int)e_last->send_message_data
		int counter; // call RecomputeFactorData() after every 10*A->arity calls to SendMessage()
	};

	void RecomputeFactorData(Energy::NonSingletonFactor* A);
};

} // namespace srmpLib
#endif
