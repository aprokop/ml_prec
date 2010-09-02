#include "solvers.h"

#include <cstdlib>
#include <cmath>

void generate_x0(Vector& x) {
#if 0
    srandom(time(NULL));
#else
    // we don't want it to change from run to run
    srandom(3);
#endif
    for (uint i = 0; i < x.size(); i++) {
	x[i] = 20.*(random() - 0.5*RAND_MAX)/RAND_MAX + 100;
	// x[i] = 1 - (i&1)*2;
	// x[i] = 1;
    }
}

double calculate_norm(const Vector& r, const CSRMatrix& A, const PrecBase& B, NormType norm_type) {
    uint n = r.size();

    switch (norm_type) {
	case NORM_L2	: return dnrm2(r);
	case NORM_A	:
	case NORM_B_1	: {
	    Vector r1(n);
	    if (norm_type == NORM_A)
		multiply(A, r, r1);
	    else {
		Vector& r_ = const_cast<Vector&>(r);
		B.solve(r_, r1);
	    }

	    double d = ddot(r1, r);
	    ASSERT(d > -1e-15, "Negative (r,r)_*: " << d);

	    return sqrt(fabs(d));
	}
    }

    THROW_EXCEPTION(" One must not be here");
}
