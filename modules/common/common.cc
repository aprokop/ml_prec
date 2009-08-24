#include "common.h"
#include "include/exception.h"

#include <cmath>


double cheb(double x, uint k) {
    ASSERT(x >= 1, "x = " << x);

    switch(k) {
	case 0: return 1;
	case 1: return x;
	case 2: return 2*x*x - 1;
	case 3: return x*(4*x*x - 3);
	default: return cosh(k*acosh(x));
    }
}
