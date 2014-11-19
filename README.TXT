SRMP, version 1.0
Vladimir Kolmogorov, vnk@ist.ac.at
http://pub.ist.ac.at/~vnk/software.html

====================================================================================================

Software for minimizing functions f(x1,...xn) of discrete variables represented as a sum of low-order terms (factors).

Implements the following algorithms:

1. SRMP - Sequential Tree-Reweighted Message Passing

	V. Kolmogorov.
	"Reweighted message passing revisited."
	arXiv report 1309.5655, September 2013.

2. CMP - Convex Message Passing (the version with uniform weights; for details see the tech report above).

	T. Hazan and A. Shashua. 
	"Norm-product belief propagation: Primal-dual message-passing for approximate inference."
	IEEE Trans. on Information Theory, 56(12):6294–6316, 2010.

3. MPLP - Max Product Linear Programming

	A. Globerson and T. Jaakkola. 
	"Fixing max-product: Convergent message passing algorithms for MAP LP-relaxations"
	NIPS, 2007.


	D. Sontag, A. Globerson and T. Jaakkola. 
	"Introduction to dual decomposition for inference."
	In S. Sra, S. Nowozin, and S. Wright, editors, Optimization for Machine Learning. MIT Press, 2011.

====================================================================================================

USAGE:

1. c++ interface - see SRMP.h.

2. Via text file. Use UAI format described in 
http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php
Example is given in grid4x4.UAI.LG . Note: it's a maximization problem!
(When reading the file, all costs are multiplied by -1).

Compile the code and run "./srmp" to see all options.

====================================================================================================

Important: for high-order factors it is possible to use different relaxations.
E.g. if there is a ternary factor ijk={i,j,k}, a pairwise factor ij={i,j} and unary factors i, j, k then
the following choices are possible:

Option 1 (weak relaxation): add edges ijk->i, ijk->j, ijk->k, ij->i, ij->k
Option 2 (stronger relaxation): add additionally edge ijk->ij; this will enforce consistency between these factors.

See "SPECIFYING RELAXATION" in SRMP.h for details.


