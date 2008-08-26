#include "vector.h"
#include "include/logger.h"

#include "config/config.h"
#ifdef HAVE_LIBBLAS
#include "include/blas.h"
#endif

#include <cmath>
#include <iostream>

DEFINE_LOGGER("Vector");

const Vector& Vector::operator=(const Vector& v) {
    // check for self-assignment
    if (this == &v)
	return *this;
    const int n = v.size();
    data.resize(v.size(), 0.);

    memcpy(&data[0], &v[0], n*sizeof(double));

    return *this;
}

Vector::Vector(const Vector& v) {
    *this = v;
}

const Vector& Vector::operator+=(const Vector& v) THROW {
    daxpy(1, v, *this);
    return *this;
}

const Vector& Vector::operator-=(const Vector& v) THROW {
    daxpy(-1, v, *this);
    return *this;
}

const Vector& Vector::operator*=(double f) {
    dscal(f, *this);
    return *this;
}

const Vector& Vector::operator/=(double f) THROW {
    ASSERT(f, "Division by zero");
    dscal(1./f, *this);
    return *this;
}

double Vector::norm_2() const {
    return dnrm2(*this);
}

// ===========================  BLAS PROCEDURES WRAPPERS  ===========================
#ifdef HAVE_LIBBLAS
void daxpy(double alpha, const Vector& x, Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    FORTRAN(daxpy)(&n, &alpha, &x[0], (int[]){1}, &y[0], (int[]){1});
}

void dscal(double alpha, Vector& x) {
    int n = x.size();

    FORTRAN(dscal)(&n, &alpha, &x[0], (int[]){1});
}

double ddot(const Vector& x, const Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    return FORTRAN(ddot)(&n, &x[0], (int[]){1}, &y[0], (int[]){1});
}

double dnrm2(const Vector& x) {
    int n = x.size();

    return FORTRAN(dnrm2)(&n, &x[0], (int[]){1});
}
#else
void daxpy(double alpha, const Vector& x, Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    for (int i = 0; i < n; i++)
	y[i] += alpha*x[i];
}

void dscal(double alpha, Vector& x) {
    int n = x.size();

    for (int i = 0; i < n; i++)
	x[i] *= alpha;
}

double ddot(const Vector& x, const Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    double s = 0.;
    for (int i = 0; i < n; i++)
	s += x[i] * y[i];
    return s;
}

double dnrm2(const Vector& x) {
    return sqrt(ddot(*this,*this));
}
#endif
