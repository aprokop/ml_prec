#include "vector.h"
#include "include/logger.h"

#include "config/config.h"
#ifdef HAVE_LIBBLAS
#include "include/blas.h"
#endif

#include <cmath>
#include <iostream>

DEFINE_LOGGER("SVector");

SVector::SVector(const Vector& v) {
    const int n = v.size();
    data.resize(v.size(), 0.);

    memcpy(&data[0], &v[0], n*sizeof(double));
}

const SVector& SVector::operator=(const SVector& v) {
    // check for self-assignment
    if (this == &v)
	return *this;
    const int n = v.size();
    data.resize(v.size(), 0.);

    memcpy(&data[0], &v[0], n*sizeof(double));

    return *this;
}

const SVector& SVector::operator+=(const SVector& v) THROW {
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

const SVector& SVector::operator-=(const SVector& v) THROW {
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

const SVector& SVector::operator*=(double f) {
    const int n = data.size();

#ifdef HAVE_LIBBLAS
    FORTRAN(dscal)(&n, &f, &data[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	data[i] *= f;
#endif

    return *this;
}

const SVector& SVector::operator/=(double f) THROW {
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
