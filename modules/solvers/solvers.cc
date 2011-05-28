#include "solvers.h"
#include "include/logger.h"

#include <cstdlib>
#include <cmath>

DEFINE_LOGGER("Solvers");

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

/* Check if the stop criteria uses too small epsilon */
const double MIN_EPS = 1e-14;
void check_and_replace_eps(double init_norm, double& eps) {
    if (init_norm*eps < MIN_EPS) {
	LOG_INFO("==========================    WARNING    ==========================");
	LOG_INFO("Bad stop criteria: ||r|| < " << init_norm*eps);
	eps = MIN_EPS/init_norm;
	LOG_INFO("New stop criteria: ||r|| < " << MIN_EPS << ".     eps = " << eps);
	LOG_INFO("==========================    WARNING    ==========================");
    }
}

std::ostream& operator<<(std::ostream& os, const SolverStats& stats) {
    os << std::fixed << std::setprecision(3);
    os << std::endl;
    os << "Residual time       :  avg = " << stats.t_resid << "\t total = " << stats.t_resid_total << std::endl;
    os << "Prec inversion time :  avg = " << stats.t_prec  << "\t total = " << stats.t_prec_total << std::endl;
    if (stats.t_const)
	os << "Time of (possible) construction: " << stats.t_const << std::endl;
    os << "Number of iter      : " << stats.niter << std::endl;
    os << "Time of solution    : " << stats.t_sol << std::endl;

    return os;
}
