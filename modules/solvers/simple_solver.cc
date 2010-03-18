#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("SimpleSolver");

void SimpleSolver(const CSRMatrix& A, const Vector& b, PrecBase& B, Vector& x,
		  double eps, bool silent) THROW {
    double  gtime = pclock();
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    uint   n = b.size();
    Vector r(n), z(n);
    double norm, init_norm;

    double  mult = 0,  inv = 0,  cstr = 0, delta;
    int	   nmult = 0, ninv = 0;

    generate_x0(x);
    residual(A, b, x, r);

    norm = init_norm = dnrm2(r);

    int niter = 0;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#endif

    while (norm/init_norm > eps) {
	if (!silent)
	    LOG_DEBUG("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	delta = pclock();
	B.solve(r, z);
	if (niter) {
	    inv += pclock() - delta;
	    ninv++;
	} else {
	    cstr = pclock() - delta;
	}

	daxpy(1., z, x);

	delta = pclock();
	residual(A, b, x, r);
	mult += pclock() - delta;
	nmult++;

	norm = dnrm2(r);

	niter++;
    }
    gtime = pclock() - gtime;

    if (!silent) {
	LLL_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	LLL_DEBUG(std::fixed << std::setprecision(3) << "Residual:       avg = " << mult/nmult << "\t total = " << mult);
	if (ninv) {
	    LLL_DEBUG(std::fixed << std::setprecision(3) << "Prec inversion: avg = " << inv/ninv << "\t total = " << inv);
	    double cstr_pos = cstr - inv/ninv;
	    if (cstr_pos > 1e-1) {
		LLL_DEBUG(std::fixed << std::setprecision(3) <<
			  "Time of (possible) construction: " << cstr_pos);
		LLL_DEBUG(std::fixed << std::setprecision(3) <<
			  "Time of (possible) solution    : " << gtime - cstr_pos);
	    }
	}
    } else {
	LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    }
}
