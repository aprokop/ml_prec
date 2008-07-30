#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <cmath>
#include <iostream>
#include <fstream>
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

#endif // __TOOLS_H__
