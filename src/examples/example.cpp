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
#include <ctype.h>
#include <math.h>
#include <srmp/SRMP.h>
#include <srmp/FactorTypes/PottsType.h>

#ifdef _MSC_VER
#pragma warning(disable: 4996) /* Disable deprecation */
#endif

using namespace srmpLib;

PottsFactorType* potts = NULL;

bool TAKE_LOG = false;

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

// Implements something similar to fgets(), but without knowing the length of strings.
// The size of internally allocated memory is at most twice the size of the longest string in the file.

class ReadLinesFromFile
{
public:
	ReadLinesFromFile(const char* filename);
	~ReadLinesFromFile();

	bool exists(); // returns true if constructor managed to open the file
	char* NextLine();
	int CurrentLine();

private:
	FILE* fp;
	char* buf;
	int size; // buf[0..size-1] contains data from the file
	int size_max; // allocated size of buf is (size_max+1)
	int size_prev; // size of previously returned line (including trailing 0)
	int line_num;
};

inline ReadLinesFromFile::ReadLinesFromFile(const char* filename)
{
	fp = fopen(filename, "rb");
	if (!fp) return;
	size_max = 2;
	buf = (char*) malloc(size_max+1);
	size = fread(buf, 1, size_max, fp);
	size_prev = 0;
	line_num = 0;
}
inline ReadLinesFromFile::~ReadLinesFromFile()
{
	if (fp) fclose(fp);
	if (buf) free(buf);
}
inline bool ReadLinesFromFile::exists()
{
	return (fp) ? true : false;
}
inline char* ReadLinesFromFile::NextLine()
{
	size -= size_prev;
	memmove(buf, buf+size_prev, size);
	size += fread(buf+size, 1, size_max-size, fp);

	if (size == 0) return NULL;

	size_prev = 0;
	while ( 1 )
	{
		for ( ; size_prev<size; size_prev++)
		{
			if (buf[size_prev] == 0 || buf[size_prev] == '\n')
			{
				buf[size_prev ++] = 0;
				line_num ++;
				return buf;
			}
		}

		if (size < size_max)
		{
			buf[size] = 0;
			line_num ++;
			return buf;
		}

		size_max *= 2;
		buf = (char*) realloc(buf, size_max+1);
		size += fread(buf+size, 1, size_max-size, fp);
	}
}
inline int ReadLinesFromFile::CurrentLine()
{
	return line_num;
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
template <typename REAL> class ReadNumbers
{
public:
	ReadNumbers();
	~ReadNumbers();
	int ReadLine(char* line); // returns # of numbers read from line
	REAL GetNumber(int i) { return arr[i]; } // can be called after ReadLine()
	REAL* GetNumbersArray() { return arr; } // can be called after ReadLine()

private:
	REAL* arr;
	int arr_size;
};

template <typename REAL> inline ReadNumbers<REAL>::ReadNumbers() : arr_size(1)
{
	arr = (REAL*) malloc(arr_size*sizeof(REAL));
}
template <typename REAL> inline ReadNumbers<REAL>::~ReadNumbers()
{
	free(arr);
}
template <typename REAL> inline int ReadNumbers<REAL>::ReadLine(char* line)
{
	REAL i;
	int num = 0;
	const char* sscanf_str = ( ((REAL)0.5) == ((REAL)0) ) ? "%d" : "%lf";
	while ( 1 )
	{
		while (line[0] && isspace(line[0])) line ++;
		if (!line[0]) return num;
		char* line_next = line+1;
		while (line_next[0] && !isspace(line_next[0])) line_next ++;
		char tmp = line_next[0];
		line_next[0] = 0;
		int sscanf_read = sscanf(line, sscanf_str, &i);
		line_next[0] = tmp;
		if (sscanf_read != 1) return num;
		line = line_next;

		num ++;
		if (num > arr_size)
		{
			arr_size *= 2;
			arr = (REAL*) realloc(arr, arr_size*sizeof(REAL));
		}
		arr[num-1] = i;
	}
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

Energy* ReadUAI(const char* filename)
{
	ReadLinesFromFile F(filename);
	if (!F.exists()) { printf("Can't open %s\n", filename); exit(1); }

	Energy* g = NULL;
	int phase = 0; // 0: reading # nodes
	               // 1: reading # states
	               // 2: reading # terms
	               // 3: reading term descriptions
	               // 4: reading term costs
	               // 5: finished
	int node_num = 0, t = 0, term_num = 0;
	struct Term
	{
		int arity;
		int* nodes;
	};
	Term* terms = NULL;
	Buffer buf(1024);

	ReadNumbers<int> read_ints;
	ReadNumbers<double> read_doubles;

	char* line;
	while ( (line=F.NextLine()) != NULL )
	{
		int i, num = read_ints.ReadLine(line);
		if (num == 0) continue;

		if (phase == 0)
		{
			if (num != 1) { printf("Error in line %d: # nodes expected\n", F.CurrentLine()); exit(1); }
			node_num = read_ints.GetNumber(0);
			if (node_num < 1) { printf("Error in line %d: # of nodes should be positive\n", F.CurrentLine()); exit(1); }
			g = new Energy(node_num);
			phase ++;
		}
		else if (phase == 1)
		{
			if (num != node_num) { printf("Error in line %d: # of nodes do not match\n", F.CurrentLine()); exit(1); }
			for (i=0; i<node_num; i++) g->AddNode(read_ints.GetNumber(i));
			phase ++;
		}
		else if (phase == 2)
		{
			if (num != 1) { printf("Error in line %d: # terms expected\n", F.CurrentLine()); exit(1); }
			term_num = read_ints.GetNumber(0);
			if (term_num < 0) { printf("Error in line %d: # of terms should be non-negative\n", F.CurrentLine()); exit(1); }
			t = 0;
			terms = new Term[term_num];
			phase ++;
		}
		else if (phase == 3)
		{
			num --;
			if (num != read_ints.GetNumber(0)) { printf("Error in line %d: incorrect term entry\n", F.CurrentLine()); exit(1); }
			if (num < 1) { printf("Error in line %d: arity must be positive\n", F.CurrentLine()); exit(1); }
			terms[t].arity = num;
			terms[t].nodes = (int*)buf.Alloc(num*sizeof(int));
			for (i=0; i<num; i++)
			{
				int q = read_ints.GetNumber(i+1);
				if (q<0 || q>=node_num)  { printf("Error in line %d: incorrect index\n", F.CurrentLine()); exit(1); }
				terms[t].nodes[i] = q;
			}
			t ++;
			if (t == term_num) { phase ++; t = 0; }
		}
		else if (phase == 4)
		{
			if (num != 1) { printf("Error in line %d: # costs expected\n", F.CurrentLine()); exit(1); }
			num = read_ints.GetNumber(0);
			line = F.NextLine();
			if (!line) { printf("Error in line %d: missing cost table\n", F.CurrentLine()); exit(1); }
			int double_num = read_doubles.ReadLine(line);
			double* arr = read_doubles.GetNumbersArray();
			if (TAKE_LOG)
			{
				for (i=0; i<double_num; i++) arr[i] = log(arr[i]);
			}
			for (i=0; i<double_num; i++) arr[i] = -arr[i]; // we are minimizing, not maximizing!
			if (num == -1)
			{
				if (terms[t].arity != 2 || g->GetK(terms[t].nodes[0]) != g->GetK(terms[t].nodes[1]))
				{ printf("Error in line %d: not a potts term\n", F.CurrentLine()); exit(1); }
				if ( 1 )
				{
					if (!potts) potts = new PottsFactorType;
					g->AddFactor(2, terms[t].nodes, arr, potts);
				}
				else
				{
					int x, y, i=0, K = g->GetK(terms[t].nodes[0]);
					double lambda = arr[0];
					double* V = new double[K*K];
					for (x=0; x<K; x++)
					for (y=0; y<K; y++)
					{
						V[i ++] = (x == y) ? 0 : lambda;
					}
					g->AddPairwiseFactor(terms[t].nodes[0], terms[t].nodes[1], V);
				}
			}
			else
			{
				int K = 1;
				for (i=0; i<terms[t].arity; i++) K *= g->GetK(terms[t].nodes[i]);
				if (double_num != K) { printf("Error in line %d: # states doesn't match\n", F.CurrentLine()); exit(1); }
				g->AddFactor(terms[t].arity, terms[t].nodes, arr);
			}
			t ++;
			if (t == term_num) phase ++;
		}
	}

	if (phase < 5)
	{
		if (g) delete g;
		g = NULL;
	}
	return g;

}


void AddTriplets(Energy* g, const char* filename)
{
	ReadLinesFromFile F(filename);
	if (!F.exists()) { printf("Can't open %s\n", filename); exit(1); }

	char* line;
	int num = 0;
	while ( (line=F.NextLine()) != NULL )
	{
		int r, i[3];
		for (r=0; r<3; r++)
		{
			while (line[0] && !isdigit(line[0])) line ++;
			if (!line[0]) break;
			char* line_next = line+1;
			while (isdigit(line_next[0])) line_next ++;
			bool stop = (line_next[0] == 0);
			line_next[0] = 0;
			i[r] = atoi(line);
			if (stop) { r ++; break; }
			line = line_next+1;
		}
		if (r==3)
		{
			//printf("Adding triplet %d,%d,%d\n", i[0], i[1], i[2]);
			g->AddTriplet(i[0], i[1], i[2]);
			num ++;
		}
	}
	printf("%d triplets from %s added\n", num, filename);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


















void ShowUsage(char* argv0, bool full_info=false)
{
	printf("Usage: %s <options> <filename.UAI.LG>\n", argv0);
	printf("general options:\n");
	printf("  -SRMP                // SRMP method\n");
	printf("  -MPLP                // MPLP method\n");
	printf("  -MPLP_BW             // MPLP with a backward pass\n");
	printf("  -CMP                 // CMP method\n");
	printf("  -iter=<n>            // maximum number of iters\n");
	printf("  -time_max=<secs>     // time limit\n");
	printf("  -eps=<v>             // stop if the increase of the lower bound during one iter is smaller than eps\n");
	if (full_info)
	{
		printf("  -t=<filename>        // tighten relaxation using triplets of nodes in the file (format: lines with \"id id id\")\n");
	}
	if (full_info)
	{
		printf("  -sort=[-1|0|1|...]   // order of processing factors. See SRMP.h for details\n");
	}
	printf("  -BLP                 // BLP relaxation (with singleton separators)\n");
	printf("  -FULL=[0|1|2]        // Full relaxation (see SetFullEdges() in SRMP.h for details)\n");
	printf("  -FULLDUAL=[...]      // Full relaxation, by constructing dual graph (see SetFullEdgesDual() in SRMP.h for details)\n");
	if (full_info)
	{
		printf("SRMP options:\n");
		printf("  -TRWS=<a>            // a\\in[0,1]. a=1 corresponds to TRW-S (for pairwise energies)\n");
		printf("\n");
		printf("  -takelog             // take logarithm after reading values from the file\n");
		printf("  -verbose=[0|1]       // 0: be silent\n");
		printf("  -save=<filename>     // save final result (primal solution) to file\n");
		printf("  -saverep=<filename>  // save final repamaterization to file\n");
		printf("  -stats               // print statistics about factors and edges\n");
		printf("\n");
		printf("To use multiple calls to Solve() with different options (e.g. first MPLP then SRMP), insert \"+\" between options\n");
	}
	else
	{
		printf("  -h                   // help: print more options\n");
	}
	exit(1);
}


struct Run
{
	Run() { tighten_filename = NULL; }
	~Run() {}

	Energy::Options options;
	char* tighten_filename;
};

int main(int argc, char* argv[])
{
	if (0)
	{
		Energy* g = ReadUAI("../grid4x4.UAI.LG");
		
		g->SetFullEdges(1);
		
		Energy::Options options;
		options.method = Energy::Options::MPLP;
		options.iter_max = 100;
		g->Solve(options);
		exit(1);
	}

	int i;

	char* filename = NULL;
	char* save_filename = NULL;
	char* save_rep_filename = NULL;
	bool BLP_relaxation = false;
	bool FULL_relaxation = false;
	bool FULLDUAL_relaxation = false;
	int FULL_relaxation_flag = 0;
	int FULLDUAL_relaxation_flag = 0;
	bool print_stats = false;
	
	Run* current = new Run;

	Block<Run*> run_list(10);
	Run** run_ptr = run_list.New();
	*run_ptr = current;

	for (i=1; i<argc; i++)
	{
		if (argv[i][0] == '+' && !argv[i][1])
		{
			current = new Run;
			run_ptr = run_list.New();
			*run_ptr = current;
			continue;
		}
		if (argv[i][0] != '-')
		{
			if (filename) { printf("Error: filename can be specified only once\n"); ShowUsage(argv[0]); }
			filename = argv[i];
			continue;
		}
		if (!strcmp("h", &argv[i][1]))
		{
			ShowUsage(argv[0], true);
		}
		if (!strncmp("save=", &argv[i][1], 5))
		{
			if (save_filename) { printf("Error: save filename can be specified only once\n"); ShowUsage(argv[0]); }
			save_filename = &argv[i][1+5];
			continue;
		}
		if (!strncmp("saverep=", &argv[i][1], 8))
		{
			if (save_filename) { printf("Error: save filename can be specified only once\n"); ShowUsage(argv[0]); }
			save_rep_filename = &argv[i][1+8];
			continue;
		}
		if (!strncmp("t=", &argv[i][1], 2))
		{
			if (current->tighten_filename) { printf("Error: tightening filename can be specified only once\n"); ShowUsage(argv[0]); }
			current->tighten_filename = &argv[i][1+2];
			continue;
		}
		if (!strncmp("iter=", &argv[i][1], 5))
		{
			current->options.iter_max = atoi(&argv[i][1+5]);
			continue;
		}
		if (!strncmp("time=", &argv[i][1], 5))
		{
			current->options.time_max = atof(&argv[i][1+5]);
			continue;
		}
		if (!strncmp("eps=", &argv[i][1], 4))
		{
			current->options.eps = atof(&argv[i][1+4]);
			continue;
		}
		if (!strcmp("SRMP", &argv[i][1]))
		{
			current->options.method = Energy::Options::SRMP;
			continue;
		}
		if (!strcmp("CMP", &argv[i][1]))
		{
			current->options.method = Energy::Options::CMP;
			continue;
		}
		if (!strcmp("MPLP", &argv[i][1]))
		{
			current->options.method = Energy::Options::MPLP;
			continue;
		}
		if (!strncmp("TRWS=", &argv[i][1], 5))
		{
			current->options.TRWS_weighting = atof(&argv[i][1+5]);
			continue;
		}
		if (!strcmp("BLP", &argv[i][1]))
		{
			BLP_relaxation = true;
			continue;
		}
		if (!strcmp("FULL", &argv[i][1]))
		{
			FULL_relaxation = true;
			continue;
		}
		if (!strncmp("FULL=", &argv[i][1], 5))
		{
			FULL_relaxation = true;
			FULL_relaxation_flag = atoi(&argv[i][1+5]);
			continue;
		}
		if (!strcmp("FULLDUAL", &argv[i][1]))
		{
			FULLDUAL_relaxation = true;
			continue;
		}
		if (!strncmp("FULLDUAL=", &argv[i][1], 9))
		{
			FULLDUAL_relaxation = true;
			FULLDUAL_relaxation_flag = atoi(&argv[i][1+9]);
			continue;
		}
		if (!strncmp("sort=", &argv[i][1], 5))
		{
			current->options.sort_flag = atoi(&argv[i][1+5]);
			continue;
		}
		if (!strncmp("verbose=", &argv[i][1], 8))
		{
			current->options.verbose = (atoi(&argv[i][1+8])) ? true : false;
			continue;
		}
		if (!strcmp("takelog", &argv[i][1]))
		{
			TAKE_LOG = true;
			continue;
		}
		if (!strcmp("stats", &argv[i][1]))
		{
			print_stats = true;
			continue;
		}
		ShowUsage(argv[0]);
	}

	if (!filename) ShowUsage(argv[0]);

	bool start = true;
	Energy* g = ReadUAI(filename);

	if (BLP_relaxation) g->SetMinimalEdges();
	if (FULL_relaxation) g->SetFullEdges(FULL_relaxation_flag);
	if (FULLDUAL_relaxation) g->SetFullEdgesDual(FULLDUAL_relaxation_flag);

	for (run_ptr=run_list.ScanFirst(); run_ptr; run_ptr=run_list.ScanNext())
	{
		current = *run_ptr;
		if (current->tighten_filename)
		{
			if (start) g->SetFullEdges();
			AddTriplets(g, current->tighten_filename);
		}
		g->Solve(current->options);
		delete current;

		start = false;
	}

	if (print_stats)
	{
		g->PrintStats();
	}

	if (save_filename)
	{
		FILE* fp = fopen(save_filename, "w");
		if (!fp) { printf("Can't open %s for writing\n", save_filename); exit(1); }
		for (i=0; i<g->GetNodeNum(); i++) fprintf(fp, "%d ", g->GetSolution(i));
		fclose(fp);
	}

	if (save_rep_filename)
	{
		g->SaveUAI(save_rep_filename, false, true);
	}

	delete g;

	return 0;
}


