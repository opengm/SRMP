/* PottsType.h */
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



/* Example usage: minimize f(x,y,z) = 2*x + x*y + 5*P(y,z)
   where x,y,z \in {0,1,2} and P(.,.) is the Potts term: P(a,b)=0 if a==b, and P(a,b)=1 otherwise.

#include <stdio.h>
#include "SRMP.h"
#include "PottsType.h"

int main(int argc, char* argv[])
{
	int i, a, b, node_num = 3;
	int f[2];
	double lambda, D[3], V[3*3];

	Energy* g = new Energy(node_num);
	PottsFactorType* potts = new PottsFactorType;

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

	// add pairwise term 5*P(y,z)
	lambda = 5;
	f[0] = 1; f[1] = 2;
	g->AddFactor(2, f, &lambda, potts);

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
	delete potts;

	return 0;
}

*/


#ifndef SRMP_NLAGNKLAJSNFRKAJSGKAS
#define SRMP_NLAGNKLAJSNFRKAJSGKAS

#include "../SRMP.h"

namespace srmpLib {

struct PottsFactorType : Energy::FactorType
{
	PottsFactorType();
	~PottsFactorType();

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
};

} // namespace srmpLib
#endif
