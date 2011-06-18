#include "solvers.h"
#include "project/config.h"
#include "modules/prec/misc/misc.h"
#include "include/logger.h"
#include "include/time.h"

DEFINE_LOGGER("ChebSolver");

void ChebSolver(const CSRMatrix& A, double lmin, double lmax, const Vector& b, const PrecBase& B, Vector& x,
		double eps, NormType norm_type, bool silent) THROW {
    double  gtime = pclock();
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    uint   n = b.size();
    Vector r(n), u0(n), u1(n);
    double alpha, beta;
    double norm, init_norm;
    double eta = (lmax + lmin) / (lmax - lmin);

    double  mult = 0,  inv = 0,  cstr = 0, delta;
    int	   nmult = 0, ninv = 0;

    generate_x0(x);
    residual(A, b, x, r);
    norm = init_norm = calculate_norm(r, A, B, norm_type);

    int niter = 1;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#else
    check_and_replace_eps(init_norm, eps);
#endif

    // ===============    STEP 1    ===============
    LOG_INFO("#" << niter-1 << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    alpha = 2/(lmax + lmin);

    u1 = x;

    // x = u0 + alpha*B.solve(b - A*u0);
    B.solve(r, x);
    dscal(alpha, x);
    daxpy(1, u1, x);
    residual(A, b, x, r);
    norm = dnrm2(r);
    niter++;

    // ===============    STEPS 2+    ===============
    while (norm/init_norm > eps) {
	if (!silent)
	    LOG_DEBUG("#" << niter-1 << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	// hack to avoid copying and allocating new memory
	u0.swap(x);
	u1.swap(u0);
	// end hack

	alpha = 4/(lmax - lmin) * cheb(eta, niter-1)/cheb(eta, niter);
	beta  = cheb(eta, niter-2) / cheb(eta, niter);

	// x = u1 + alpha*B.solve(b - A*u1) + beta*(u1 - u0);
	delta = pclock();
	B.solve(r, x);
	inv += pclock() - delta;
	ninv++;

	for (uint k = 0; k < n; k++)
	    x[k] = u1[k] + alpha*x[k] + beta*(u1[k] - u0[k]);

	delta = pclock();
	residual(A, b, x, r);
	mult += pclock() - delta;
	nmult++;
	norm = calculate_norm(r, A, B, norm_type);

	niter++;
    }
    gtime = pclock() - gtime;

    niter--;

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
