#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("SimpleSolver");

void SimpleSolver(const CSRMatrix& A, const Vector& b, PrecBase& B, Vector& x,
		  double eps, bool silent) THROW {
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    uint   n = b.size();
    Vector r(n), z(n);
    double norm, init_norm;

    clock_t  mult = 0,  inv = 0,  cstr = 0, delta;
    int	    nmult = 0, ninv = 0;

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

	delta = clock();
	B.solve(r, z);
	if (niter) {
	    inv += clock() - delta;
	    ninv++;
	} else {
	    cstr = clock() - delta;
	}

	daxpy(1., z, x);

	delta = clock();
	residual(A, b, x, r);
	mult += clock() - delta;
	nmult++;

	norm = dnrm2(r);

	niter++;
    }

    if (!silent) {
	double mult_ = double(mult) / CLOCKS_PER_SEC;
	double inv_  = double(inv)  / CLOCKS_PER_SEC;

	LLL_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	LLL_DEBUG(std::fixed << std::setprecision(3) << "Residual:       avg = " << mult_/nmult << "\t total = " << mult_);
	if (ninv) {
	    LLL_DEBUG(std::fixed << std::setprecision(3) << "Prec inversion: avg = " << inv_/ninv << "\t total = " << inv_);
	    LLL_DEBUG(std::fixed << std::setprecision(3) <<
		      "Time of (possible construction) [time of first inversion - avg] = " <<
		      double(cstr)/CLOCKS_PER_SEC - inv_/ninv);
	}
    } else {
	LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    }
}
