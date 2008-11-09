#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("PCG");

Vector PCGSolver(const CSRMatrix& A, const Vector& b, PrecBase& B, double eps) THROW {
    ASSERT(A.rows() == A.cols() && A.cols() == b.size(), "Wrong dimesions: " << 
	   "A:" << A.rows() << " x " << A.cols() << ", b: " << b.size());

    uint   n = b.size();
    Vector r(n), x(n), z(n);
    Vector p(n), Ap(n);
    double alpha, beta;
    double app, rz0, rz1;
    double norm, init_norm;

    clock_t  mult = 0,  inv = 0,  cstr = 0, delta;
    int	    nmult = 0, ninv = 0;

    LOG_DEBUG("Generating initial approximation");
    generate_x0(x);
    residual(A, b, x, r);
    norm = init_norm = r.norm_2();

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
	LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	delta = clock();
	multiply(A, p, Ap);
	mult += clock() - delta;
	nmult++;

	app = ddot(Ap, p);
	if (is_equal(app, 0.)) {
	    LOG_WARN("ddot(Ap, p) = 0");
	    break;
	}
	alpha = rz0 / app;

	daxpy(alpha, p, x);
	daxpy(-alpha, Ap, r);

	norm = r.norm_2();

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
