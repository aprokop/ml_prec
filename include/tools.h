#ifndef __TOOLS_H__
#define __TOOLS_H__

#include "include/define.h"
#include "include/uvector.h"

#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <sys/types.h>
#include <signal.h>

#define FORTRAN(f) f ## _

static const double epsilon = std::numeric_limits<double>::epsilon();

inline bool is_equal(double x, double y) {
    return (fabs(x-y) <= epsilon) ? true : false;
}

inline bool is_not_equal(double x, double y) {
    return (fabs(x-y) > epsilon) ? true : false;
}

template<typename T>
inline T random(T i0, T i1) {
    assert(i0 < i1);
    return i0 + (i1-i0)*(double(random())/RAND_MAX);
}

template<typename T>
const std::vector<T> std_vector(uint n, T t1, ...) {
    std::vector<T> v(n);

    v[0] = t1;
    va_list ap;
    va_start(ap, t1);
    for (uint i = 0; i < n-1; i++)
	v[i+1] = va_arg(ap, T);
    va_end(ap);

    return v;
}

template<typename T>
const uvector<T> new_uvector(uint n, T t1, ...) {
    uvector<T> v(n);

    v[0] = t1;
    va_list ap;
    va_start(ap, t1);
    for (uint i = 0; i < n-1; i++)
	v[i+1] = va_arg(ap, T);
    va_end(ap);

    return v;
}


template<typename T>
inline T pow2(const T& x) {
    return x*x;
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
