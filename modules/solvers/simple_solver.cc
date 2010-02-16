#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("SimpleSolver");

Vector SimpleSolver(const CSRMatrix& A, const Vector& b, PrecBase& B, double eps) THROW {
    ASSERT_SIZE(b.size(), A.size());

    uint   n = b.size();
    Vector r(n), x(n), z(n);
    double norm, init_norm;

    clock_t  mult = 0,  inv = 0,  cstr = 0, delta;
    int	    nmult = 0, ninv = 0;

    LOG_DEBUG("Generating initial approximation");
    generate_x0(x);
    residual(A, b, x, r);

    norm = init_norm = dnrm2(r);

    int niter = 0;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#endif

    while (norm/init_norm > eps) {
	LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

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

    LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    std::cout << "#" << niter << ": relative -> " << std::scientific << norm/init_norm << std::fixed << std::endl;

#if 1
    double out = double(mult)/CLOCKS_PER_SEC;
    LLL_DEBUG(std::fixed << std::setprecision(3) << "Residual:       avg = " << out/nmult << "\t total = " << out);
    out = double(inv)/CLOCKS_PER_SEC;
    if (ninv) {
	LLL_DEBUG(std::fixed << std::setprecision(3) << "Prec inversion: avg = " << out/ninv << "\t total = " << out);
	LLL_DEBUG(std::fixed << std::setprecision(3) <<
		  "Time of (possible construction) [time of first inversion - avg] = " <<
		  double(cstr - inv/ninv)/CLOCKS_PER_SEC);
    }
#endif

    return x;
}
