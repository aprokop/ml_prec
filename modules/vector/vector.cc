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
    const int n = data.size();
    ASSERT(v.size() == uint(n), "Adding vector with different dimension: " << v.size() << " (" << n << ")");

#ifdef HAVE_LIBBLAS
    FORTRAN(daxpy)(&n, (double[]){1}, &v[0], (int[]){1}, &data[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	data[i] += v[i];
#endif

    return *this;
}

const Vector& Vector::operator-=(const Vector& v) THROW {
    const int n = data.size();
    ASSERT(v.size() == uint(n), "Subtracting vector with different dimension: " << v.size() << " (" << n << ")");

#ifdef HAVE_LIBBLAS
    FORTRAN(daxpy)(&n, (double[]){-1}, &v[0], (int[]){1}, &data[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	data[i] -= v[i];
#endif

    return *this;
}

const Vector& Vector::operator*=(double f) {
    const int n = data.size();

#ifdef HAVE_LIBBLAS
    FORTRAN(dscal)(&n, &f, &data[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	data[i] *= f;
#endif

    return *this;
}

const Vector& Vector::operator/=(double f) THROW {
    ASSERT(f, "Division by zero");
    const int n = data.size();

#ifdef HAVE_LIBBLAS
    double d = 1./f;
    FORTRAN(dscal)(&n, &d, &data[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	data[i] /= f;
#endif

    return *this;
}

double Vector::norm_2() const {
#ifdef HAVE_LIBBLAS
    int n = data.size();
    return FORTRAN(dnrm2)(&n, &data[0], (int[]){1});
#else
    return sqrt(scalar_product(*this,*this));
#endif
}

double scalar_product(const Vector& v1, const Vector& v2) THROW {
    const uint n = v1.size();
    ASSERT(v2.size() == n, "Trying to find scalar product of vectors with different dimensions");

#ifdef HAVE_LIBBLAS
    int _n = n;
    return FORTRAN(ddot)(&_n, &v1[0], (int[]){1}, &v2[0], (int[]){1});
#else
    double s = 0.;
    for (int i = 0; i < n; i++)
	s += v1[i] * v2[i];
    return s;
#endif
}


