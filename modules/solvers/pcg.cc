#include "solvers.h"
#include "include/logger.h"
#include "include/time.h"

#ifdef HAVE_LIBBLAS
#include "include/blas.h"
#endif

DEFINE_LOGGER("PCG");

Vector PCG(const MatrixInterface& A, const Vector& b, PrecBase& B, double eps) THROW {
    ASSERT(A.rows() == A.cols() && A.cols() == b.size(), "Wrong dimesions: " << 
	   "A:" << A.rows() << " x " << A.cols() << ", b: " << b.size());

    int    n = b.size();
    Vector r(n), x(n), z(n);
    Vector p(n), Ap(n);
    double alpha, beta;
    double rz0, rz1;
    double norm, init_norm;

    clock_t  mult = 0,  inv = 0,  cstr = 0, delta;
    int	    nmult = 0, ninv = 0;

#if 1
    LOG_DEBUG("Generating random initial approximation");
    srandom(time(NULL)); 
    for (int i = 0; i < n; i++) {
	x[i] = 20.*(random() - 0.5*RAND_MAX)/RAND_MAX + 100;
	// x[i] = 400.*(random() - 0.5*RAND_MAX)/RAND_MAX + 300;
	// x[i] = 1 - (i&1)*2;
	// x[i] = 100*sin(3*i);
	// x[i] = 1;
	// x[i] = 100;
    }
    delta = clock();
    r = b - A*x;
    mult += clock() - delta;
    nmult++;
#else
    // We take x0 = 0
    r = b;
#endif

    delta = clock();
    // LOG_DEBUG("r = " << std::scientific << r);
    B.solve(r, z);
    cstr = clock() - delta;

    memcpy(&p[0], &z[0], n*sizeof(double));
    // LOG_DEBUG("p = " << std::scientific << p);
    rz0 = scalar_product(r, z);
    norm = init_norm = r.norm_2();
    // exit(0);

    int niter = 1;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#endif

    while (norm/init_norm > eps) {
	LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	delta = clock();
	// LOG_DEBUG("p = " << std::scientific << p);
	multiply(A, p, Ap);
	// LOG_DEBUG("Ap = " << std::scientific << Ap);
	mult += clock() - delta;
	nmult++;

	alpha = rz0 / scalar_product(Ap, p);

#ifdef HAVE_LIBBLAS
	FORTRAN(daxpy)(&n, &alpha, &p[0], (int[]){1}, &x[0], (int[]){1});
	double nalpha = -alpha;
	FORTRAN(daxpy)(&n, &nalpha, &Ap[0], (int[]){1}, &r[0], (int[]){1});
#else
	x += alpha*p;
	r -= alpha*Ap;
#endif
	// LOG_DEBUG("x = " << std::scientific << x);
	// LOG_DEBUG("r = " << std::scientific << r);

	norm = r.norm_2();

	delta = clock();
	B.solve(r, z);
	inv += clock() - delta;
	ninv++;

	rz1 = scalar_product(r, z);
	beta = rz1 / rz0;

#ifdef HAVE_LIBBLAS
	FORTRAN(dscal)(&n, &beta, &p[0], (int[]){1});
	FORTRAN(daxpy)(&n, (double[]){1}, &z[0], (int[]){1}, &p[0], (int[]){1});
#else
	p = z + beta*p;
#endif

	rz0 = rz1;
	niter++;
    }

    LOG_INFO("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    std::cout << "#" << niter << ": relative -> " << std::scientific << norm/init_norm << std::fixed << std::endl;

#if 1
#undef LOG_DEBUG
#define LOG_DEBUG(v)	 std::cout << "DEBUG : " << __func__ << " : " << v << std::endl 

    double out = double(mult)/CLOCKS_PER_SEC;
    LOG_DEBUG(std::fixed << std::setprecision(3) << "Residual:       avg = " << out/nmult << "\t total = " << out);
    out = double(inv)/CLOCKS_PER_SEC;
    LOG_DEBUG(std::fixed << std::setprecision(3) << "Prec inversion: avg = " << out/ninv << "\t total = " << out);
    LOG_DEBUG(std::fixed << std::setprecision(3) << 
	      "Time of (possible construction) [time of first inversion - avg] = " << 
	      double(cstr - inv/ninv)/CLOCKS_PER_SEC);
#endif

    return x;
}
