/* SharedPairwiseType.h */
/* Copyright Vladimir Kolmogorov vnk@ist.ac.at */
/* Created May 2013 */

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

/* 
   Example usage: minimize f(x,y,z) = 2*x + 4*V(x,y) + 5*V(y,z)
   where x,y,z \in {0,1,2} and V(.,.) is the following function: V(a,b)=a*b.

   In this example the same term V is reused several times (with different positive coefficients).
   The advantage of using SharedPairwiseFactorType is memory efficiency:
   we need to allocate term V only once.

#include <stdio.h>
#include "SRMP.h"
#include "SharedPairwiseType.h"

int main(int argc, char* argv[])
{
	int i, a, b, node_num = 3;
	int f[2];
	double lambda, D[3], V[3*3];

	Energy* g = new Energy(node_num);

	for (a=0; a<3; a++)
	for (b=0; b<3; b++)
	{
		V[3*a + b] = a*b;
	}
	SharedPairwiseFactorType* shared_pairwise = new SharedPairwiseFactorType(3, 3, V);

	for (i=0; i<node_num; i++)
	{
		g->AddNode(3); // add node with 3 labels
	}

	// add unary term 2*x
	D[0] = 0; D[1] = 2; D[2] = 4;
	g->AddUnaryFactor(0, D);

	// add pairwise term 4*V(x,y)
	f[0] = 0; f[1] = 1;
	lambda = 4;
	g->AddFactor(2, f, &lambda, shared_pairwise);

	// add pairwise term 5*V(y,z)
	f[0] = 1; f[1] = 2;
	lambda = 5;
	g->AddFactor(2, f, &lambda, shared_pairwise);

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
	delete shared_pairwise;

	return 0;
}

*/

#ifndef SRMP_KAJSFKJAFKJASGAKJSFAAGARARAGS
#define SRMP_KAJSFKJAFKJASGAKJSFAAGARARAGS

#include "../SRMP.h"

namespace srmpLib {

struct SharedPairwiseFactorType : Energy::FactorType
{
	SharedPairwiseFactorType(int K1, int K2, double* costs);
	~SharedPairwiseFactorType();

	void InitFactor(Energy::NonSingletonFactor* A, double* user_data);
	double GetCost(Energy::NonSingletonFactor* A);

	void ComputePartialReparameterization(Energy::NonSingletonFactor* A, double* rep);

	void InitEdge(Energy::Edge* e);

	double SendMessage(Energy::Edge* e);
	void SendRestrictedMessage(Energy::Edge* e);

	double SendMPLPMessages(Energy::NonSingletonFactor* A, bool set_solution);

	bool PrepareFactor(Energy::NonSingletonFactor* A);

	/////////////////////////////////////////////////////
private:
	int K1, K2;
	double* costs;
	Buffer buf;
	ReusableBuffer rbuf;
	ReusableBuffer rbuf2;
};

} // namespace srmpLib
#endif
