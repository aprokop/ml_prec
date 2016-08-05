#ifndef __TOOLS_H__
#define __TOOLS_H__

#include "include/define.h"
#include "include/uvector.h"

#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <signal.h>

static const double epsilon = std::numeric_limits<double>::epsilon();

/*
 * Because of floating numbers are approximations to real numbers, we do not
 * use strict comparison for double. Instead, we allow to numbers to be declared
 * the same if their difference is small enough
 */
inline bool is_equal(double x, double y) {
    return (fabs(x-y) <= epsilon) ? true : false;
}

inline bool is_not_equal(double x, double y) {
    return (fabs(x-y) > epsilon) ? true : false;
}

/*
 * Create a random number of type T in the interval [i0, i1].
 * NOTE: random(-1,1) would return an integer; to get double call random(-1.0,1.0)
 */
template<typename T>
inline T random(T i0, T i1) {
    assert(i0 < i1);
    return i0 + (i1-i0)*(double(random())/RAND_MAX);
}

/*
 * Create std::vector from a list of scalars of type T.
 * The first argument is the number of variables in the list.
 * Typical use:
 *	std::vector<double> x = new_vector(4, 1.0, 0.0, -2, 4);
 */
template<typename T>
std::vector<T> new_vector(uint n, T t1, ...) {
    std::vector<T> v(n);

    v[0] = t1;
    va_list ap;
    va_start(ap, t1);
    for (uint i = 0; i < n-1; i++)
        v[i+1] = va_arg(ap, T);
    va_end(ap);

    return v;
}

/*
 * Create std::vector from a string.
 * Typical use:
 *	std::vector<double> x = new_vector<double>("1.0 0.0 -2 4");
 */
template<typename T>
std::vector<T> new_vector(const std::string& str) {
    std::vector<T> v;

    std::istringstream is(str);
    std::copy(std::istream_iterator<T>(is), std::istream_iterator<T>(), std::back_inserter(v));

    return v;
}

/* The same as std_vector but for uvector */
template<typename T>
uvector<T> new_uvector(uint n, T t1, ...) {
    uvector<T> v(n);

    v[0] = t1;
    va_list ap;
    va_start(ap, t1);
    for (uint i = 0; i < n-1; i++)
        v[i+1] = va_arg(ap, T);
    va_end(ap);

    return v;
}

/* Compute the square of x */
template<typename T>
inline T pow2(const T& x) {
    return x*x;
}

/* Check if a number is NaN */
template<typename T>
inline bool is_nan(const T& x) {
    return (((x) != (x)) || (((x) < (x))));
}

#if 0
/*
 * Macros to be used when wants to measure memory output.
 * Throws SIGUSR1 signal to be caught by the runner
 */
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
