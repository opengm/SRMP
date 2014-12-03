/*
    PQ.h - implements pairing heaps priority queue (with multipass 'delete min')

    Copyright 2008 Vladimir Kolmogorov (vnk@ist.ac.at)

    This software can be used for research purposes only. Commercial use is prohibited.
    Public redistribution of the code or its derivatives is prohibited.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef SRMP_HFKSJHFKJHARBABDAKFAF
#define SRMP_HFKSJHFKJHARBABDAKFAF

// exactly one flag must be defined
//#define SRMP_PQ_MULTIPASS
#define SRMP_PQ_INTERLEAVED_MULTIPASS

#include <string.h>

namespace srmpLib {

template <typename REAL> class PriorityQueue
{
public:
	struct Item
	{
		REAL	slack;

		Item*	parentPQ;
		union Children
		{
			struct Direction
			{
				Item*	leftPQ;
				Item*	rightPQ;
			} direction;
			REAL	y_saved; // used in repairs
		} children;
	};
	static void* AllocateBuf();
	static void DeallocateBuf(void* buf);

	static void ResetItem(Item* i);
	static bool isReset(Item* i);

	//////////////////////////////////////////////////////////

	void Reset();
	void Add(Item* i);
#define SRMP_Remove(i, buf) _Remove(i)
	void _Remove(Item* i);
	void Decrease(Item* i_old, Item* i_new/*, void* buf*/);
	Item* GetMin();

	//////////////////////////////////////////////////////////

	void Update(REAL delta);
	void Merge(PriorityQueue<REAL>& dest);

	// traversing items in the order they are stored (irrespective of slack).
	// The caller must go through all items, no other member functions can be called during the scan.
	Item* GetAndResetFirst();
	Item* GetAndResetNext();

	Item* GetFirst();
	Item* GetNext(Item* i);

	//////////////////////////////////////////////////////////

private:
	struct Buf
	{
	};
	Item*	rootPQ;
	void RemoveRoot();
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline void* PriorityQueue<REAL>::AllocateBuf()
{
	return NULL;
}

template <typename REAL> inline void PriorityQueue<REAL>::DeallocateBuf(void* _buf)
{
}

template <typename REAL> inline void PriorityQueue<REAL>::ResetItem(Item* i) 
{ 
	i->parentPQ = NULL;
}

template <typename REAL> inline bool PriorityQueue<REAL>::isReset(Item* i) 
{ 
	return (i->parentPQ == NULL);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline void PriorityQueue<REAL>::Reset() 
{ 
	rootPQ = NULL; 
}

/*
template <typename REAL> inline void PriorityQueue<REAL>::RemoveRoot()
{
	Item* r = rootPQ;
	PriorityQueue<REAL> pq;
	pq.rootPQ = rootPQ;
	rootPQ = NULL;
	Item* i;
	for (i=pq.GetAndResetFirst(); i; i=pq.GetAndResetNext())
	{
		if (i != r) Add(i);
	}
	r->parentPQ = NULL;
}
*/

// sets i = merge(i, j). Ignores parentPQ and rightPQ for i and j.
#define SRMP_MERGE_PQ(i, j)\
	{\
		if (i->slack <= j->slack)\
		{\
			j->children.direction.rightPQ = i->children.direction.leftPQ;\
			if (j->children.direction.rightPQ) j->children.direction.rightPQ->parentPQ = j;\
			j->parentPQ = i;\
			i->children.direction.leftPQ = j;\
		}\
		else\
		{\
			i->children.direction.rightPQ = j->children.direction.leftPQ;\
			if (i->children.direction.rightPQ) i->children.direction.rightPQ->parentPQ = i;\
			i->parentPQ = j;\
			j->children.direction.leftPQ = i;\
			i = j;\
		}\
	}

template <typename REAL> inline void PriorityQueue<REAL>::RemoveRoot()
{
	Item* i = rootPQ->children.direction.leftPQ;
	rootPQ->parentPQ = NULL;
	if (i)
	{
#ifdef SRMP_PQ_MULTIPASS
		while ( i->children.direction.rightPQ )
		{
			Item** prev_ptr = &rootPQ;
			while ( 1 )
			{
				if (i->children.direction.rightPQ)
				{
					Item* j = i->children.direction.rightPQ;
					Item* next = j->children.direction.rightPQ;
					SRMP_MERGE_PQ(i, j);
					*prev_ptr = i;
					if (!next) { i->children.direction.rightPQ = NULL; break; }
					prev_ptr = &i->children.direction.rightPQ;
					i = next;
				}
				else
				{
					*prev_ptr = i;
					i->children.direction.rightPQ = NULL;
					break;
				}
			}
			i = rootPQ;
		}
#endif

#ifdef SRMP_PQ_INTERLEAVED_MULTIPASS
		while ( i->children.direction.rightPQ )
		{
			Item* prev = NULL;
			while ( i )
			{
				Item* next;
				if (i->children.direction.rightPQ)
				{
					Item* j = i->children.direction.rightPQ;
					next = j->children.direction.rightPQ;
					SRMP_MERGE_PQ(i, j);
				}
				else next = NULL;
				i->children.direction.rightPQ = prev;
				prev = i;
				i = next;
			}
			i = prev;
		}
#endif
		i->parentPQ = i;
	}
	rootPQ = i;
}

template <typename REAL> inline void PriorityQueue<REAL>::Add(Item* i)
{
	if (!rootPQ)
	{
		rootPQ = i;
		i->parentPQ = i;
		i->children.direction.leftPQ = i->children.direction.rightPQ = NULL;
	}
	else if (i->slack <= rootPQ->slack)
	{
		rootPQ->parentPQ = i;
		i->children.direction.leftPQ = rootPQ;
		i->children.direction.rightPQ = NULL;
		rootPQ = i;
		i->parentPQ = i;
	}
	else
	{
		i->children.direction.leftPQ = NULL;
		i->children.direction.rightPQ = rootPQ->children.direction.leftPQ;
		if (i->children.direction.rightPQ) i->children.direction.rightPQ->parentPQ = i;
		rootPQ->children.direction.leftPQ = i;
		i->parentPQ = rootPQ;
	}
}


template <typename REAL> inline void PriorityQueue<REAL>::_Remove(Item* i)
{
	Item* p = i->parentPQ;
	if (p == i) RemoveRoot();
	else
	{
		if (i->children.direction.rightPQ) i->children.direction.rightPQ->parentPQ = p;
		if (p->children.direction.leftPQ == i) p->children.direction.leftPQ  = i->children.direction.rightPQ;
		else                p->children.direction.rightPQ = i->children.direction.rightPQ;
		if (i->children.direction.leftPQ)
		{
			i->parentPQ = i;
			i->children.direction.rightPQ = NULL;
			PriorityQueue<REAL> pq;
			pq.rootPQ = i;
			pq.RemoveRoot();
			pq.Merge(*this);
		}
		else i->parentPQ = NULL;
	}
}

template <typename REAL> inline void PriorityQueue<REAL>::Decrease(Item* i_old, Item* i_new/*, void* _buf*/)
{
	if (i_old->parentPQ == i_old)
	{
		if (i_old != i_new)
		{
			rootPQ = i_new;
			i_new->parentPQ = i_new;
			i_new->children.direction.leftPQ = i_old->children.direction.leftPQ;
			i_new->children.direction.rightPQ = NULL;
			if (i_new->children.direction.leftPQ) i_new->children.direction.leftPQ->parentPQ = i_new;
			i_old->parentPQ = NULL;
		}
	}
	else
	{
		SRMP_Remove(i_old, _buf);
		Add(i_new);
	}
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetMin()
{
	return rootPQ;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



template <typename REAL> inline void PriorityQueue<REAL>::Merge(PriorityQueue<REAL>& dest)
{
	if (!rootPQ) return;
	if (!dest.rootPQ) dest.rootPQ = rootPQ;
	else
	{
		if (rootPQ->slack < dest.rootPQ->slack)
		{
			Item* j = rootPQ; rootPQ = dest.rootPQ; dest.rootPQ = j;
		}
		rootPQ->children.direction.rightPQ = dest.rootPQ->children.direction.leftPQ;
		if (rootPQ->children.direction.rightPQ) rootPQ->children.direction.rightPQ->parentPQ = rootPQ;
		rootPQ->parentPQ = dest.rootPQ;
		dest.rootPQ->children.direction.leftPQ = rootPQ;
	}
	rootPQ = NULL;
}



template <typename REAL> inline void PriorityQueue<REAL>::Update(REAL delta)
{
	if (!rootPQ) return;

	Item* i = rootPQ;
	while (i->children.direction.leftPQ) i = i->children.direction.leftPQ;

	while ( 1 )
	{
		i->slack += delta;

		if (i->children.direction.rightPQ)
		{
			i = i->children.direction.rightPQ;
			while (i->children.direction.leftPQ) i = i->children.direction.leftPQ;
		}
		else
		{
			while ( 1 )
			{
				Item* j = i;
				i = i->parentPQ;
				if (i == j) return;
				if (i->children.direction.leftPQ == j) break;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetAndResetFirst()
{
	if (!rootPQ) return NULL;
	return GetAndResetNext();
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetAndResetNext()
{
	if (!rootPQ) return NULL;
	Item* result = rootPQ;
	result->parentPQ = NULL;
	Item* i = rootPQ->children.direction.leftPQ;
	if (!i) rootPQ = result->children.direction.rightPQ;
	else
	{
		rootPQ = i;
		while (i->children.direction.rightPQ) i = i->children.direction.rightPQ;
		i->children.direction.rightPQ = result->children.direction.rightPQ;
	}
	return result;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetFirst()
{
	if (!rootPQ) return NULL;
	Item* i = rootPQ;
	while (i->children.direction.leftPQ) i = i->children.direction.leftPQ;
	return i;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetNext(Item* i)
{
	if (i->children.direction.rightPQ)
	{
		i = i->children.direction.rightPQ;
		while (i->children.direction.leftPQ) i = i->children.direction.leftPQ;
		return i;
	}
	while ( 1 )
	{
		Item* j = i;
		i = i->parentPQ;
		if (i == j) return NULL;
		if (i->children.direction.leftPQ == j) return i;
	}
}

} // namespace srmpLib

#endif
