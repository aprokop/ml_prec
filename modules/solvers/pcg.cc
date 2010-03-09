#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("PCG");

void PCGSolver(const CSRMatrix& A, const Vector& b, PrecBase& B, Vector& x, double eps, bool silent) THROW {
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    uint   n = b.size();
    Vector r(n), z(n);
    Vector p(n), Ap(n);
    double alpha, beta;
    double app, rz0, rz1;
    double norm, init_norm;

    clock_t  mult = 0,  inv = 0,  cstr = 0, delta;
    int	    nmult = 0, ninv = 0;

    generate_x0(x);
    residual(A, b, x, r);
    norm = init_norm = dnrm2(r);

#if 0
    // w-Jacobi relaxation (w = 0.5)
    double w = 0.5;
    for (uint i = 0; i < n; i++)
	x[i] += r[i] / (w*A.get(i,i));
    residual(A, b, x, r);
    for (uint i = 0; i < n; i++)
	x[i] += r[i] / (w*A.get(i,i));
    residual(A, b, x, r);
#endif

    delta = clock();
    B.solve(r, z);
    cstr = clock() - delta;

    p = z;
    rz0 = ddot(r, z);

    int niter = 0;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#endif

    while (norm/init_norm > eps) {
	if (!silent)
	    LOG_DEBUG("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	delta = clock();
	multiply(A, p, Ap);
	mult += clock() - delta;
	nmult++;

	app = ddot(Ap, p);
	if (app <= 0) {
	    LOG_ERROR("<Ap, p> <= 0");
	    THROW_EXCEPTION("<Ap, p> = " << app << ": |Ap| = " << dnrm2(Ap) << ", |p| = " << dnrm2(p));
	}
	alpha = rz0 / app;

	daxpy(alpha, p, x);
	daxpy(-alpha, Ap, r);

	norm = dnrm2(r);

	delta = clock();
	B.solve(r, z);
	inv += clock() - delta;
	ninv++;

	rz1 = ddot(r, z);
	beta = rz1 / rz0;

	dscal(beta, p);
	daxpy(1, z, p);

	rz0 = rz1;
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
