/* SRMP.h */
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
   Software for minimizing functions  f(x1,...,xn) of discrete variables x1, ..., xn
   represented as a sum of low-order terms (or "factors").

   Hopefully, the usage should be self-explanatory. For concrete examples see FactorTypes/PottsType.h and FactorTypes/SharedPairwiseType.h.
*/

#ifndef FNAIJDGAKSJDBAJFGBASFJAS
#define FNAIJDGAKSJDBAJFGBASFJAS

#include <string.h>
#include <assert.h>
#include "Algs/block.h"
#include "Algs/timer.h"

//#define VNK_DEBUG


// flags for Energy::AddFactor(...)
const unsigned FLAG_DO_NOT_COPY_INTO_INTERNAL_MEMORY = 0x0001; // used for PairwiseFactorType, GeneralFactorType and PatternFactorType

#define ENERGY_INFTY (1e100)

class Energy
{
public:
	typedef int NodeId;
	typedef void* FactorId;
	struct FactorType;

	Energy(int node_num_max);
	~Energy();

	/////////////////////////////////////////////////////////////////////////
	//                         CREATING INSTANCE                           //
	/////////////////////////////////////////////////////////////////////////

	// first call returns 0, second returns 1, and so on. K = number of labels.
	//
	// If costs!=NULL then it points to the array of size 'K'; 
	// the pointer is just copied, no new memory for costs is allocated
	NodeId AddNode(int K, double* costs=NULL);

	// 'costs' is an array of size GetK(i); it is added to the array in internal memory
	FactorId AddUnaryFactor(NodeId i, double* costs);
	// 'costs' is an array of size GetK(i)*GetK(j), with cost(xi,xj)=costs[GetK(j)*xi + xj].
	// It is copied into internal memory
	FactorId AddPairwiseFactor(NodeId i, NodeId j, double* costs);

	// Custom factors and factors of arbitrary arity.
	//
	// nodes is an array of size 'arity', where arity>=1
	// costs is an array of size nodes[0].K * ... * nodes[arity-1].K.
	// The last node corresponds to the "least significant bit" (as in AddPairwiseFactor(), with nodes[0]=i and nodes[1]=j)
	//
	// costs can also be NULL (= zero cost array)
	//
	// Use 'type' for custom factor types; see examples in PottsType.h, PatternType.h and SharedPairwiseType.h.
	// The meaning of pointer 'costs' may be different for custom types (and the same for 'flags')
	FactorId AddFactor(int arity, NodeId* node_indexes, double* costs, FactorType* type = NULL, unsigned flags = 0);

	int GetNodeNum() { return node_num; }
	int GetK(NodeId i) { return nodes[i].K; }


	/////////////////////////////////////////////////////////////////////////
	//                         CALLING SOLVER                              //
	/////////////////////////////////////////////////////////////////////////

	struct Options
	{
		Options()
			: method(SRMP),
			  iter_max(10000),
			  time_max(20*60), // 20 minutes
			  eps(1e-8),
			  compute_solution_period(10),
			  print_times(true),
			  sort_flag(0),
			  verbose(true),
			  TRWS_weighting(1.0)
		{}

		enum
		{
			SRMP,
			MPLP,
			MPLP_BW, // MPLP with a backward pass
			CMP
		} method;

		int iter_max;
		double time_max; // in seconds
		double eps; // stop if the increase of the lower during one iteration is less than 'eps'
		int compute_solution_period; // extract solution after every 'compute_solution_period' iterations
		bool print_times;
		int sort_flag; // -1: process factors in the order they were given (except that nodes in SRMP and CMP are always traversed first)
		               //  0: sort factors according to the given node ordering
		               //  1: use an automatic greedy technique for sorting nodes, then sort factors accordingly
		               // >1: user random permutation of nodes, with 'sort_flag' as the seed
		bool verbose;

		// SRMP options
		double TRWS_weighting; // in [0,1], 1 corresponds to TRW-S (for pairwise energies)
	};

	double Solve(Options& options);
	int GetSolution(NodeId i) { return nodes[i].solution_best; } // can be called after Solve()


	/////////////////////////////////////////////////////////////////////////
	//                         SPECIFYING RELAXATION                       //
	/////////////////////////////////////////////////////////////////////////

	// The relaxation is specified by 'edges', i.e. pairs of factors (A,B) where B is a strict subset of A.
	//
	// Functions SetMinimalEdges() and SetFullEdges() can be called only if no edges have been added yet.
	// If no edges have been added before calling Solve() then Solve() calls SetFullEdges().
	//
	// If you add more factors after calling Solve() then you'll need to add edges for them explicitly.

	void SetMinimalEdges(); // for each factor A={i1,...,ik} add edges (A,{i1}), ..., (A,{ik})

	// method=0: add all possible pairs (A,B) with B \subset A, with the following exception:
	//           if there exists factor C with B\subset C \subset A then don't add (A,B)
	// method=1: move all costs to outer factors (converting them to general types first, if they are not already of these types).
	//           Then run method=0 and remove unnecesary edges (i.e. those that do not affect the relaxation).
	//           Note, all edges outgoing from outer factors will be kepts.
	// method=2: similar to method=1, but all edges {A->B, B->C} are replaced with {A->B, A->C} (so this results in a two-layer graph).
	//
	// Note, method=1 and method=2 merge duplicate factors while method=0 does not. For this reason the relaxation may be tighther.
	// (If there are no duplicate factors then the resulting relaxation should be the same in all three cases).
	void SetFullEdges(int method=0);

	// method=3: run method=2 and then create a new Energy instance with unary and pairwise terms in which nodes correspond
	// to outer factors of the original energy, and pairwise terms with {0,+\infty} costs enforce consistency between them.
	// 'sort_flag' has the same meaning as in Options::sort_flag.
	// NOTE: IF THIS FUNCTION IS CALLED THEN AFTERWARDS ONLY   Solve() and GetSolution()   ARE ALLOWED TO BE CALLED!
	void SetFullEdgesDual(int sort_flag);


	// Specify edges manually. B must be a subset of A
	void AddRelaxationEdge(FactorId A, FactorId B);
	static bool isSuperset(FactorId A, FactorId B); // for sanity, may want to check first that B is indeed a subset of A

	// adds a factor on {i,j,k} with zero costs and edges to factors {i,j}, {j,k}, {i,k}.
	// If such pairwise factors do not exist they are added first (with edges to singleton factors, if add_all_edges is true)
	void AddTriplet(NodeId i, NodeId j, NodeId k, bool add_all_edges=false);
	
	
	
	/////////////////////////////////////////////////////////////////////////
	//                             MISCELLANEOUS                           //
	/////////////////////////////////////////////////////////////////////////
	// save in UAI.LG format. 
	// if sort_factors is false then store factors in the order that they were added
	void SaveUAI(char* filename, bool sort_factors=false, bool save_reparameterization=false);
	void PrintStats(); // print # of factors and edges
	
	
	
	
	
	
	
	
	
	
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

// all structures are public, to allow adding custom factor types

	struct Factor;
	struct Node;
	struct NonSingletonFactor;
	struct Edge;

	struct Factor
	{
		int arity;
		int K;
		double* data; // the meaning of 'data' depends on the type
		double* rep; // this var can be used arbitrarily by different methods. E.g. MPLP uses it to store reparameterization for separators.

		// flags
		unsigned int is_unsorted : 1; // if is_unsorted==1 then nodes+arity contains a sorted list of nodes (if arity > 1)
		unsigned int is_removed : 1;
		unsigned int tmp1 : 1;
		unsigned int tmp2 : 1;
		//unsigned int is_degenerate_node

		// 3 variables below are valid only for separators
		unsigned int compute_bound : 1; // = 1 if this factor is last touched during the backward pass
		unsigned int weight_forward : 21;
		unsigned int weight_backward : 21;

		unsigned int user1 : 1; // can be used by implementations of FactorType

		Edge* first_in;
	};
	struct Node : Factor 
	{
		int solution, solution_best;
	};
	struct NonSingletonFactor : Factor
	{
		Node** nodes; // of size arity
		FactorType* type;
		Edge* first_out;

		// theta is an array of size 'K'. Minimize theta over labelings consistent with Node::solution
		// (for nodes i with i->solution >= 0). Set Node::solution for unlabeled nodes appropriately.
		// _buf must be an array of size 4*arity.
		void ComputeRestrictedMinimum(double* theta, int* _buf);
		// similar to the above, but instead of setting Node::solution set thetaB (where factor B is a subset of factor 'this')
		void ComputeRestrictedMinimum(double* theta, Factor* B, double* thetaB, int* _buf);
	};
	struct Edge
	{
		NonSingletonFactor* A; // B must be a subset of A
		Factor* B;
		Edge* next_out;
		Edge* next_in;

		unsigned int is_fw : 1; // = 0 is B is the last factor among A+children of A
		unsigned int is_bw : 1; // = 0 if B is the first factor among A+children of A

		unsigned int compute_bound : 1; // = 1 if factor A is last touched during the backward pass
		unsigned int weight_forward : 21;   // in the current implementation
		unsigned int weight_backward : 21;  // these weights are 0 or 1

		unsigned int user1 : 1; // can be used by implementations of FactorType
		unsigned int user2 : 1; // can be used by implementations of FactorType
		void* send_message_data; // can be used by implementations of FactorType

		double* m; // message, of size B->K.   m=NULL means that A->type->InitEdge() has not been called yet.
	};
	struct FactorType
	{
// if TEST_FACTOR_TYPES is defined then outputs of ComputePartialRepameterization(), SendMessage(), SendRestrictedMessage(), SendMPLPMessages()
// will be compared against true outputs (where costs are obtained via GetCost()). 
// Exception: the set_solution option in SendMPLPMessages() is not tested.
//
//#define TEST_FACTOR_TYPES 

		// called whenever a new factor is added
		virtual void InitFactor(NonSingletonFactor* A, double* user_data, unsigned flags) = 0;

		// return the cost of the labeling stored in Node::solution
		virtual double GetCost(NonSingletonFactor* A) = 0;

		// set   theta[xA] = costs[xA] - \sum_{(A,C)} m_{AC}[xC]
		// where the sum is over edges outgoing from A, m_{AC} is the message from A to C, and xC is the restriction of xA to C.
		// (A 'full' reparameterization would also include messages along incoming edges, but in this function they are ignored.)
		//
		// This function would not be if incoming edges are disallowed in PrepareFactor().
		virtual void ComputePartialReparameterization(NonSingletonFactor* A, double* theta) = 0;

		// called whenever a new outgoing edge is added
		virtual void InitEdge(Edge* e) = 0; // there holds e->A->type == this

		// 1. for edge e=(A,B) set
		//   e->m[xB] = min_{xA~xB} { A->costs[xA] + \sum_{(D,A)} m_{DA}[xA] - \sum_{(A,C)!=(A,B)} m_{AC}[xC] }
		// where xA~xB means that labelings xA and xB are consistent (i.e. the restriction of xA to B is xB)
		// 2. compute delta = min_xB e->m[xB]
		//    normalize e->m[xB] -= delta
		// 3. return delta
		//
		// If the code calls SendMessage(e1) and SendMessage(e2), then between these two calls only e1->m may change;
		// other messages are guaranteed to stay the same. This could be exploited for speeding up computations;
		// see example in PatternType.h
		virtual double SendMessage(Edge* e) = 0; // there holds e->A->type == this

		// similar to SendMessage, but the minimum is computed only over labelings xA consistent with currently labeled nodes
		// (which is indicated by i->solution for Node i: i->labeled<0 means unlabeled, i->labeled>=0 means labeled).
		// An arbitrary constant can be added to all entries of e->m.
		// This function is used for extracting a primal solution. It is called only if e->B contains unlabeled nodes
		// and e->A contains labeled nodes.
		virtual void SendRestrictedMessage(Edge* e) = 0; // there holds e->A->type == this

		// 1. set  theta[xA] = A->costs[xA] + \sum_{(D,A)} m_{DA}[xA] + \sum_{(A,B)} B->rep[xB] 
		// 2. compute delta = min_xA theta[xA]
		//    normalize theta[xA] -= delta
		// 3. if set_solution is true set unlabeled nodes in A to minimize theta[xA]
		// 4. for each (A,B) set   B->rep[xB] = rho_{AB} \min_{xA~xB} theta[xA]
		//    where rho_{AB} = weight_forward_{AB} / ( weight_forward_A + \sum_{AB} weight_forward_{AB} )
		// 5. if A->rep != NULL (which happens when A has incoming edges) then set
		//      A->rep[xA] = theta[xA] - \sum_{(A,B)} B->rep[xB]
		// 6. return delta
		virtual double SendMPLPMessages(NonSingletonFactor* A, bool set_solution=false) = 0;

		// This function is called just before starting a message passing algorithm.
		// Potential use: it might inspect incoming/outgoing edges and switch to a different type, if necessary,
		// or initialize some data structure based on messages.
		//
		// If the function returns false the factor will be converted to either 'PairwiseType' or 'GeneralType'
		virtual bool PrepareFactor(NonSingletonFactor* A) = 0;
	};











private:
	int node_num, node_num_max, arity_max;
	Node* nodes;
	Block<NonSingletonFactor> factors;
	Block<Edge>* edges;
	void AllocateEdges() { if (!edges) edges = new Block<Edge>(512); }
	Buffer buf;

	static Node** GetSortedNodesPtr(NonSingletonFactor* A)
	{
		return A->nodes + A->is_unsorted*A->arity;
	}
	static Node** GetSortedNodesPtr(Factor*& A)
	{
		if (A->arity == 1) return (Node**)&A;
		else               return GetSortedNodesPtr((NonSingletonFactor*)A);
	}
	static int CompareFactors(Factor* A, Factor* B);
	static int CompareFactorsX(Factor* A, Factor* B); // use Node::solution as the node's key
	static bool isSuperset(NonSingletonFactor* A, NonSingletonFactor* B);
	static void IntersectFactors(NonSingletonFactor* A, NonSingletonFactor* B, NonSingletonFactor* C);

	struct FactorPtr
	{
		Factor* A;
		bool operator<(const FactorPtr& B) const { return (CompareFactors(A, B.A)<0); }
		bool operator>(const FactorPtr& B) const { return (CompareFactors(A, B.A)>0); }
	};
	struct FactorPtrX
	{
		Factor* A;
		bool operator<(const FactorPtrX& B) const { return (CompareFactorsX(A, B.A)<0); }
		bool operator>(const FactorPtrX& B) const { return (CompareFactorsX(A, B.A)>0); }
	};
	struct Sequence // sequence of factors
	{
		int num;
		FactorPtr* arr; // sometimes will be cast to FactorPtrX*
		int arity_max, K_max;
	};

	bool is_solution_best_initialized;
	bool is_cost_best_valid;
	double cost_best;

	// ordering
	void SortSequence(Sequence& seq, Options& options);
	void GenerateNodeOrdering(Options& options);
	int ScoreNodeOrdering(); // of ordering stored in Node::solution
	Node* _GenerateNodeOrdering(Node* seed); // write it to Node::solution. return the last node. CAN ONLY BE CALLED FROM GenerateNodeOrdering() !!!

	/////////////////////////////////////////////////////////
	FactorType* factor_type_general;
	FactorType* factor_type_pairwise;
	void InitFactor(NonSingletonFactor* A, double* user_data, FactorType* type = NULL, unsigned flags = 0);

	void GetFactorCosts(NonSingletonFactor* A, double* costs, void* _buf); // theta is an array of size A->K, _buf is array of size A->arity*sizeof(int)
	void ConvertFactorToGeneral(NonSingletonFactor* A, void* _buf); // _buf is array of size A->K*sizeof(double) + A->arity*sizeof(int)

	void SetMaximalEdges(); // performs SetFullEdges(0)
	void MoveCostsToOuterFactors();

	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////

	void InitEdges(); // allocates messages and calls FactorType::InitEdge

	// InitXXX function return the initial value of the lower bound
	double InitSRMP(Sequence& seq, Options& options);
	double SolveSRMP(Options& options);
	double InitMPLP(Sequence& seq, Options& options); // MPLP + MPLP_BW
	double SolveMPLP(Options& options);
	double InitCMP(Sequence& seq, Options& options);
	double SolveCMP(Options& options);


	void ComputeSolution(Factor* B, void* _buf); // called for separators during Solve(). _buf is of size 2*B->K*sizeof(double) + 4*B->arity*sizeof(int)

	double ComputeCost(); // of labeling stored in Node::solution
	double ComputeLowerBound(); // mainly for debugging


	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	// dual graph computations

	FactorType* factor_type_pairwise_dual;
	Energy* primal_graph;
	Energy* dual_graph;
	Sequence* dual_sequence;
	bool dual_solution_was_inconsistent;

	double SolveDualGraph(Options& options);
	double ConvertSolutionDualToPrimal();
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
/*
	// tightening - TODO
	struct Triplet
	{
		void GetNodes(Node*& i, Node*& j, Node*& k) { i=e[0]->nodes[0]; j = e[0]->nodes[1]; k = (e[1]->nodes[0]==i || e[1]->nodes[0]==j) ? e[1]->nodes[1] : e[1]->nodes[0]; }
		NonSingletonFactor* e[3]; // three edges
		double gap;
	};
	Block<Triplet> triplets;
public:
	void FindTriplets(); // finds all possible triplets, adds them to 'triplets'
	void Tighten();
*/
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	// for debugging
public:
	void PrintFactor(Factor* A) { Node** Anodes = GetSortedNodesPtr(A); int i; for (i=0; i<A->arity; i++) printf("%d", (int)(Anodes[i]-nodes)); }
	void PrintEdge(Edge* e) { PrintFactor(e->A); printf("->"); PrintFactor(e->B); }
	void Print();

#ifdef TEST_FACTOR_TYPES 
	void   TestComputePartialReparameterization(NonSingletonFactor* A, double* theta);
	double TestSendMessage(Edge* e);
	void   TestSendRestrictedMessage(Edge* e);
	double TestSendMPLPMessages(NonSingletonFactor* A, bool set_solution=false);
#endif

#ifdef VNK_DEBUG
	void AddRandomEdges(double prob);
	FactorId GetFactorId(int arity, NodeId* node_indexes); // returns the first such factor (or NULL, if such factor doesn't exist)
#endif
};

inline int Energy::CompareFactors(Factor* A, Factor* B)
{
	Node** Anodes = GetSortedNodesPtr(A);
	Node** Bnodes = GetSortedNodesPtr(B);
	int r, rA = A->arity, rB = B->arity;

	if (Anodes[0] < Bnodes[0]) return -1;
	if (Anodes[0] > Bnodes[0]) return +1;

	if (Anodes[rA-1] < Bnodes[rB-1]) return -1;
	if (Anodes[rA-1] > Bnodes[rB-1]) return +1;

	if (rA < rB) return -1;
	if (rA > rB) return +1;

	for (r=1; r<rA-1; r++)
	{
		if (Anodes[r] < Bnodes[r]) return -1;
		if (Anodes[r] > Bnodes[r]) return +1;
	}
	return 0;
}

inline int Energy::CompareFactorsX(Factor* A, Factor* B)
{
	Node** Anodes = GetSortedNodesPtr(A);
	Node** Bnodes = GetSortedNodesPtr(B);
	int r, rA = A->arity, rB = B->arity;

	if (Anodes[0]->solution < Bnodes[0]->solution) return -1;
	if (Anodes[0]->solution > Bnodes[0]->solution) return +1;

	if (Anodes[rA-1]->solution < Bnodes[rB-1]->solution) return -1;
	if (Anodes[rA-1]->solution > Bnodes[rB-1]->solution) return +1;

	if (rA < rB) return -1;
	if (rA > rB) return +1;

	for (r=1; r<rA-1; r++)
	{
		if (Anodes[r]->solution < Bnodes[r]->solution) return -1;
		if (Anodes[r]->solution > Bnodes[r]->solution) return +1;
	}
	return 0;
}







////////////////////////////////////////////////////////////////





#ifdef TEST_FACTOR_TYPES

#define COMPUTE_PARTIAL_REPARAMETERIZATION(A, theta) TestComputePartialReparameterization(A, theta) 
#define SEND_MESSAGE(e)                              TestSendMessage(e)
#define SEND_RESTRICTED_MESSAGE(e)                   TestSendRestrictedMessage(e)
#define SEND_MPLP_MESSAGES(A, set_solution)          TestSendMPLPMessages(A, set_solution)

#else

#define COMPUTE_PARTIAL_REPARAMETERIZATION(A, theta) (A)->type->ComputePartialReparameterization(A, theta) 
#define SEND_MESSAGE(e)                              (e)->A->type->SendMessage(e)
#define SEND_RESTRICTED_MESSAGE(e)                   (e)->A->type->SendRestrictedMessage(e)
#define SEND_MPLP_MESSAGES(A, set_solution)          (A)->type->SendMPLPMessages(A, set_solution)

#endif

#endif

