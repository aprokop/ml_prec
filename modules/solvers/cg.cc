#include "solvers.h"
#include "include/logger.h"

#ifdef HAVE_LIBBLAS
#include "include/blas.h"
#endif

DEFINE_LOGGER("CG");

Vector CG(const MatrixInterface& A, const Vector& b, double eps) THROW {
    ASSERT(A.rows() == A.cols() && A.cols() == b.size(), "Wrong dimesions: " << 
	   "A:" << A.rows() << " x " << A.cols() << ", b: " << b.size());

    int    n = b.size();
    Vector r(n), x(n), z(n);
    Vector p(n), Ap(n);
    double alpha, beta;
    double rz0, rz1;
    double norm, init_norm;

    p = r = b; 
    norm = init_norm = r.norm_2();
    rz0 = norm*norm;

    int niter = 1;
    while(norm > eps) {
	// if (!(niter & 0xff))
	    // LOG_DEBUG("#" << niter-1 << ": " << std::scientific << norm);
	Ap = A*p;

	alpha = rz0 / scalar_product(Ap, p);

#ifdef HAVE_LIBBLAS
	FORTRAN(daxpy)(&n, &alpha, &p[0], (int[]){1}, &x[0], (int[]){1});
	double nalpha = -alpha;
	FORTRAN(daxpy)(&n, &nalpha, &Ap[0], (int[]){1}, &r[0], (int[]){1});
#else
	x += alpha*p;
	r -= alpha*Ap;
#endif

	norm = r.norm_2();
	rz1 = norm*norm;
	beta = rz1 / rz0;

#ifdef HAVE_LIBBLAS
	FORTRAN(dscal)(&n, &beta, &p[0], (int[]){1});
	FORTRAN(daxpy)(&n, (double[]){1}, &r[0], (int[]){1}, &p[0], (int[]){1});
#else
	p = r + beta*p;
#endif

	rz0 = rz1;
	niter++;
    }
    // LOG_DEBUG("#" << niter-1 << ": " << std::scientific << norm);

    return x;
}
