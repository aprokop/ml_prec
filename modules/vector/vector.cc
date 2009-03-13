#include "vector.h"
#include "include/logger.h"
#include "include/tools.h"

#include "config/config.h"
#ifdef HAVE_BLAS
#include "include/blas.h"
#endif

#include <cmath>
#include <cstring>
#include <iostream>

DEFINE_LOGGER("Vector");

const Vector& Vector::operator=(const Vector& v) {
    // check for self-assignment
    if (this == &v)
	return *this;
    const int n = v.size();
    data.resize(n, 0.);

    if (n)
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

bool Vector::operator>=(double v) const {
    const uint n = data.size();
    for (uint i = 0; i < n; i++)
	if (data[i] < v)
	    return false;
    return true;
}

double Vector::norm_2() const {
    return dnrm2(*this);
}

bool Vector::is_nan() const {
    const uint n = data.size();
    for (uint i = 0; i < n; i++)
	if (::is_nan(data[i]))
	    return true;
    return false;
}


// ===========================  BLAS PROCEDURES WRAPPERS  ===========================
#ifdef HAVE_BLAS
void daxpy(double alpha, const Vector& x, Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    int one = 1;
    FORTRAN(daxpy)(&n, &alpha, &x[0], &one, &y[0], &one);
}

void dscal(double alpha, Vector& x) {
    int n = x.size();

    int one = 1;
    FORTRAN(dscal)(&n, &alpha, &x[0], &one);
}

double ddot(const Vector& x, const Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    int one = 1;
    return FORTRAN(ddot)(&n, &x[0], &one, &y[0], &one);
}

double dnrm2(const Vector& x) {
    int n = x.size();

    int one = 1;
    return FORTRAN(dnrm2)(&n, &x[0], &one);
}
#else
void daxpy(double alpha, const Vector& x, Vector& y) {
    int n = x.size();
    ASSERT((int)y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");

    if (!is_equal(alpha, 1.))
	for (int i = 0; i < n; i++)
	    y[i] += alpha*x[i];
    else
	for (int i = 0; i < n; i++)
	    y[i] += x[i];
}

void dscal(double alpha, Vector& x) {
    int n = x.size();

    if (is_equal(alpha, 1.))
	return;
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
    return sqrt(ddot(x,x));
}
#endif
