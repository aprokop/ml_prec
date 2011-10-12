#include "vector.h"
#include "include/blas.h"
#include "include/tools.h"

// ===========================  BLAS PROCEDURES WRAPPERS  ===========================
#ifdef HAVE_BLAS
void daxpy(double alpha, const double* x, double* y, uint n) {
    int one = 1;
    int in  = n;
    FORTRAN(daxpy)(&in, &alpha, x, &one, y, &one);
}

void dscal(double alpha, double* x, uint n) {
    int one = 1;
    int in  = n;
    FORTRAN(dscal)(&in, &alpha, x, &one);
}

double ddot(const double* x, const double* y, uint n) {
    int one = 1;
    int in  = n;
    return FORTRAN(ddot)(&in, &x[0], &one, &y[0], &one);
}

double dnrm2(const double* x, uint n) {
    int one = 1;
    int in  = n;
    return FORTRAN(dnrm2)(&in, x, &one);
}
#else
void daxpy(double alpha, const double* x, double* y, uint n) {
    if (!is_equal(alpha, 1.))
	for (uint i = 0; i < n; i++)
	    y[i] += alpha*x[i];
    else
	for (uint i = 0; i < n; i++)
	    y[i] += x[i];
}

void dscal(double alpha, double* x, uint n) {
    for (uint i = 0; i < n; i++)
	x[i] *= alpha;
}

double ddot(const double* x, const double* y, uint n) {
    double s = 0.;
    for (uint i = 0; i < n; i++)
	s += x[i] * y[i];
    return s;
}

double dnrm2(const double* x, uint n) {
    return sqrt(ddot(x,x,n));
}
#endif

void daxpy(double alpha, const Vector& x, Vector& y) {
    uint n = x.size();
    ASSERT(y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    if (n)
	daxpy(alpha, &x[0], &y[0], n);
}

void dscal(double alpha, Vector& x) {
    int n = x.size();

    if (n && !is_equal(alpha, 1.))
	dscal(alpha, &x[0], n);
}

double ddot(const Vector& x, const Vector& y) {
    uint n = x.size();
    ASSERT(y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    return ddot(&x[0], &y[0], n);
}

double dnrm2(const Vector& x) {
    int n = x.size();
    return dnrm2(&x[0], n);
}

Vector vector_product(const Vector& v1, const Vector& v2) THROW {
    ASSERT(v1.size() == 3, "Wrong v1 size: " << v1.size());
    ASSERT(v2.size() == 3, "Wrong v2 size: " << v2.size());

    Vector v(3);
    v[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return v;
}

std::ostream& operator<<(std::ostream& os, const Vector& v) {
    const uint n = v.size();
    os << "Size = " << n << std::endl;
    for (uint i = 0; i < n; i++)
	os << "   " << i << ": " << v[i] << std::endl;
    return os;
}
