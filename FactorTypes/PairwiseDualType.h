/* PairwiseDualType.h */
/* Copyright Vladimir Kolmogorov vnk@ist.ac.at */
/* Created July 2013 */

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

#ifndef LAKSJFLAKSFLAKSFJNALSFFAISF
#define LAKSJFLAKSFLAKSFJNALSFFAISF

#include "../SRMP.h"

struct PairwiseDualFactorType : Energy::FactorType
{
	PairwiseDualFactorType();
	~PairwiseDualFactorType();

	// user_data will be cast to Energy::Edge**
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

	static int ComputeRestriction(Energy::Edge* e, int k);

	Buffer buf;
	ReusableBuffer rbuf;

};


#endif
