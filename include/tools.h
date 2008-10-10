#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <sys/types.h>
#include <signal.h>

typedef unsigned int uint;

#define FORTRAN(f) f ## _

static const double epsilon = std::numeric_limits<double>::epsilon();

inline bool is_equal(double x, double y) {
    return (fabs(x-y) <= epsilon) ? true : false;
}

inline bool is_not_equal(double x, double y) {
    return (fabs(x-y) > epsilon) ? true : false;
}

inline uint random(uint i0, uint i1) {
    assert(i0 < i1);
    return i0 + uint((i1-i0)*double(random())/RAND_MAX);
}

template<typename T>
inline bool is_nan(const T& x) {
    return (((x) != (x)) || (((x) < (x))));
}

#if 0
#define LEAVE_MESSAGE(msg) { \
    pid_t ppid = getppid(); \
    std::ofstream os(".msg"); \
    if (os.good()) { \
	os << msg << std::endl; \
	kill(ppid, SIGUSR1); \
    } \
}
#else
#define LEAVE_MESSAGE(ms) 
#endif

#endif // __TOOLS_H__
