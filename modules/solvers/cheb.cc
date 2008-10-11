#include "solvers.h"
#include "modules/common/common.h"
#include "include/logger.h"
#include "include/time.h"

#include <cstdlib>

DEFINE_LOGGER("ChebSolver");

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

Vector ChebSolver(const CSRMatrix& A, double lmin, double lmax, 
		  const Vector& b, PrecBase& B, double eps) THROW {
    ASSERT(A.rows() == A.cols() && A.cols() == b.size(), "Wrong dimesions: " << 
	   "A:" << A.rows() << " x " << A.cols() << ", b: " << b.size());

    uint   n = b.size();
    Vector r(n), x(n), u0(n), u1(n);
    double alpha, beta;
    double norm, init_norm;
    double eta = (lmax + lmin) / (lmax - lmin);

    clock_t  mult = 0,  inv = 0,  cstr = 0, delta;
    int	    nmult = 0, ninv = 0;

    LOG_DEBUG("Generating initial approximation");
    generate_x0(x);
    residual(A, b, x, r);
    norm = init_norm = r.norm_2();

    int niter = 1;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#endif

    // ===============    STEP 1    ===============
    LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    alpha = 2/(lmax + lmin);

    u1 = x;

    // x = u0 + alpha*B.solve(b - A*u0);
    B.solve(r, x);
    dscal(alpha, x);
    daxpy(1, u1, x);
    residual(A, b, x, r);
    norm = r.norm_2();
    niter++;

    // ===============    STEPS 2+    ===============
    while (norm/init_norm > eps) {
	LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	// hack to avoid copying and allocating new memory
	u0.swap(x);
	u1.swap(u0);
	// end hack

	alpha = 4/(lmax - lmin) * cheb(eta, niter-1)/cheb(eta, niter);
	beta  = cheb(eta, niter-2) / cheb(eta, niter);

	// x = u1 + alpha*B.solve(b - A*u1) + beta*(u1 - u0);
	delta = clock();
	B.solve(r, x);
	inv += clock() - delta;
	ninv++;

	for (uint k = 0; k < n; k++) 
	    x[k] = u1[k] + alpha*x[k] + beta*(u1[k] - u0[k]);

	delta = clock();
	residual(A, b, x, r);
	mult += clock() - delta;
	nmult++;
	norm = r.norm_2();

	niter++;
    }

    LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    std::cout << "#" << niter << ": relative -> " << std::scientific << norm/init_norm << std::fixed << std::endl;

#if 1
#define LLL_DEBUG(v)	 LOG_DEBUG(v);std::cout << v << std::endl 
    double out = double(mult)/CLOCKS_PER_SEC;
    LLL_DEBUG(std::fixed << std::setprecision(3) << "Residual:       avg = " << out/nmult << "\t total = " << out);
    out = double(inv)/CLOCKS_PER_SEC;
    LLL_DEBUG(std::fixed << std::setprecision(3) << "Prec inversion: avg = " << out/ninv << "\t total = " << out);
    LLL_DEBUG(std::fixed << std::setprecision(3) << 
	      "Time of (possible construction) [time of first inversion - avg] = " << 
	      double(cstr - inv/ninv)/CLOCKS_PER_SEC);
#undef LLL_DEBUG
#endif

    return x;
}
