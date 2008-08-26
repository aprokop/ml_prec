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
#ifdef HAVE_LIBBLAS
    int n = data.size();
    return FORTRAN(dnrm2)(&n, &data[0], (int[]){1});
#else
    return sqrt(ddot(*this,*this));
#endif
}

void daxpy(double alpha, const Vector& x, Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

#ifdef HAVE_LIBBLAS
    FORTRAN(daxpy)(&n, &alpha, &x[0], (int[]){1}, &y[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	y[i] += alpha*x[i];
#endif
}

double ddot(const Vector& x, const Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

#ifdef HAVE_LIBBLAS
    return FORTRAN(ddot)(&n, &x[0], (int[]){1}, &y[0], (int[]){1});
#else
    double s = 0.;
    for (int i = 0; i < n; i++)
	s += x[i] * y[i];
    return s;
#endif
}

void dscal(double alpha, Vector& x) {
    int n = x.size();

#ifdef HAVE_LIBBLAS
    FORTRAN(dscal)(&n, &alpha, &x[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	x[i] *= alpha;
#endif
}
