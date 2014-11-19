/* timer.h */
/* Copyright Vladimir Kolmogorov vnk@ist.ac.at */
/* Copied from the 'Blossom V' package */



#ifndef SRMP_NJAKSDTHASKJERAXJGFBZJDLAGZ
#define SRMP_NJAKSDTHASKJERAXJGFBZJDLAGZ

// At most one of the flags
//    SRMP_PM_TIMER_MSVC
//    SRMP_PM_TIMER_CLOCK_GETTIME
//    SRMP_PM_TIMER_GETRUSAGE
//    SRMP_PM_TIMER_EXTERNAL
//    SRMP_PM_TIMER_NONE
// can be defined externally. If SRMP_PM_TIMER_EXTERNAL is defined,
// then there must exist a definition of function "double get_time()" elsewhere.

#if defined (SRMP_PM_TIMER_MSVC) || defined (SRMP_PM_TIMER_CLOCK_GETTIME) || defined (SRMP_PM_TIMER_GETRUSAGE) || defined (SRMP_PM_TIMER_EXTERNAL) || defined (SRMP_PM_TIMER_NONE)
#else
	// default option
	#ifdef _MSC_VER
		#define SRMP_PM_TIMER_MSVC
	#elif defined(__APPLE_CC__)
		#define SRMP_PM_TIMER_GETRUSAGE
	#else
		#define SRMP_PM_TIMER_CLOCK_GETTIME
	#endif
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef SRMP_PM_TIMER_MSVC

	#include <windows.h>

    namespace srmpLib {

    inline double get_time()
	{
		LARGE_INTEGER t, frequency;
		QueryPerformanceCounter(&t);
		QueryPerformanceFrequency(&frequency);
		return (double)t.QuadPart/(double)frequency.QuadPart;
	}

    } // namespace srmpLib
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef SRMP_PM_TIMER_CLOCK_GETTIME

	#include <time.h>

    namespace srmpLib {

	inline double get_time()
	{
		struct timespec t;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
		return (double)t.tv_nsec*1.00E-9 + (double)t.tv_sec;
	}

    } // namespace srmpLib
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef SRMP_PM_TIMER_GETRUSAGE

	#include <sys/resource.h>

    namespace srmpLib {

	inline double get_time()
	{
		struct rusage t;
		getrusage (RUSAGE_SELF, &t);
		return (double)t.ru_utime.tv_usec*1.00E-6 + (double)t.ru_utime.tv_sec;
	}

    } // namespace srmpLib
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef SRMP_PM_TIMER_EXTERNAL

    namespace srmpLib {

	extern double get_time();

    } // namespace srmpLib
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef SRMP_PM_TIMER_NONE

    namespace srmpLib {

	inline double get_time() { return 0; }

    } // namespace srmpLib
#endif

#endif

