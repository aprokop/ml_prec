#include "vector.h"
#include "include/logger.h"

#ifdef HAVE_LIBBLAS
#include "include/blas.h"
#endif

#include <cmath>
#include <iostream>


DEFINE_LOGGER("Vector");

Vector::Vector(uint _n) {
    n = _n;
    data.reset(new data_type(n, 0));
}

Vector::Vector(const std::vector<double>& v) {
    n = v.size();
#ifdef HAVE_LIBBLAS
    data.reset(new data_type(n));
    int _n = n;
    FORTRAN(dcopy)(&_n, &v[0], (int[]){1}, &(*data)[0], (int[]){1});
#else
    data.reset(new data_type(v));
#endif
}

Vector::Vector(const Vector& v) {
    n = v.n;
    data = v.data;
}

uint Vector::size() const {
    return n;
}

void Vector::resize(uint new_size, double v) {
    set_unique();
    n = new_size;
    data->resize(n, v);
}

double& Vector::operator[](uint index) THROW {
    check_index(index);
    set_unique();

    return (*data)[index];
}

bool Vector::operator==(const Vector& v) const THROW {
    ASSERT(n == v.n, "Comparing vectors with different dmensions");
    return (*data) == (*v.data);
}

bool Vector::operator!=(const Vector& v) const THROW {
    return !(*this == v);
}

const Vector& Vector::operator+() const {
    return *this;
}

Vector Vector::operator-() const {
    Vector v(n);
    for (uint i = 0; i < n; i++)
	(*v.data)[i] = -(*data)[i];
    return v;
}

Vector Vector::operator+(const Vector& v) const THROW {
    ASSERT(v.size() == n, "Summing up vectors with different dimensions");

    Vector w(n);
#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(dcopy)(&_n, &(*data)[0], (int[]){1}, &w[0], (int[]){1});
    FORTRAN(daxpy)(&_n, (double[]){1}, &v[0], (int[]){1}, &w[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	w[i] = (*data)[i] + v[i];
#endif
    return w;
}

Vector Vector::operator-(const Vector& v) const THROW {
    ASSERT(v.size() == n, "Summing up vectors with different dimensions");

    Vector w(n);
#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(dcopy)(&_n, &(*data)[0], (int[]){1}, &w[0], (int[]){1}); 
    FORTRAN(daxpy)(&_n, (double[]){-1}, &v[0], (int[]){1}, &w[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	w[i] = (*data)[i] - v[i];
#endif
    return w;
}

Vector Vector::operator*(double f) const {
    Vector v(n);
#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(dcopy)(&_n, &(*data)[0], (int[]){1}, &v[0], (int[]){1});
    FORTRAN(dscal)(&_n, &f, &v[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	v[i] = f*(*data)[i];
#endif
    return v;
}

Vector Vector::operator/(double f) const THROW {
    ASSERT(f, "Division by zero");

    Vector v(n);
#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(dcopy)(&_n, &(*data)[0], (int[]){1}, &v[0], (int[]){1});
    double d = 1./f;
    FORTRAN(dscal)(&_n, &d, &v[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	v[i] = (*data)[i] / f;
#endif
    return v;
}

const Vector& Vector::operator=(const Vector& v) {
    n = v.n;
    data = v.data;

    return *this;
}

void Vector::copy(const Vector& v) {
    ASSERT(n == v.n, "Trying to copy vector with different dimensions");
    memcpy(&(*data)[0], &(*v.data)[0], n*sizeof(data_type));
}

const Vector& Vector::operator+=(const Vector& v) THROW {
    ASSERT(v.n == n, "Adding vector with different dimension");
    set_unique();

#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(daxpy)(&_n, (double[]){1}, &v[0], (int[]){1}, &(*data)[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	(*data)[i] += v[i];
#endif

    return *this;
}

const Vector& Vector::operator-=(const Vector& v) THROW {
    ASSERT(v.n == n, "Subtracting vector with different dimension");
    set_unique();

#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(daxpy)(&_n, (double[]){-1}, &v[0], (int[]){1}, &(*data)[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	(*data)[i] -= v[i];
#endif

    return *this;
}

const Vector& Vector::operator*=(double f) {
    set_unique();

#ifdef HAVE_LIBBLAS
    int _n = n;
    FORTRAN(dscal)(&_n, &f, &(*data)[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	(*data)[i] *= f;
#endif

    return *this;
}

const Vector& Vector::operator/=(double f) THROW {
    ASSERT(f, "Division by zero");
    set_unique();

#ifdef HAVE_LIBBLAS
    double d = 1./f;
    int _n = n;
    FORTRAN(dscal)(&_n, &d, &(*data)[0], (int[]){1});
#else
    for (int i = 0; i < n; i++)
	(*data)[i] /= f;
#endif

    return *this;
}

double Vector::norm_1() const {
#ifdef HAVE_LIBBLAS
    int _n = n;
    return FORTRAN(dasum)(&_n, &(*data)[0], (int[]){1});
#else
    double s = 0.;
    for (int i = 0; i < n; i++)
	s += fabs((*data)[i]);
    return s;
#endif
}

double Vector::norm_2() const {
#ifdef HAVE_LIBBLAS
    int _n = n;
    return FORTRAN(dnrm2)(&_n, &(*data)[0], (int[]){1});
#else
    return sqrt(scalar_product(*this,*this));
#endif
}

double Vector::norm_inf() const {
#ifdef HAVE_LIBBLAS
    int _n = n;
    return FORTRAN(idamax)(&_n, &(*data)[0], (int[]){1});
#else 
    double s = 0.;
    for (int i = 0; i < n; i++)
	if (fabs((*data)[i]) > s) s = fabs((*data)[i]);
    return s;
#endif
}

std::ostream& operator<<(std::ostream& os, const Vector& v) {
    os << "Size = " << v.size() << std::endl;
    for (uint i = 0; i < v.size(); i++)
	os << "   " << i << ": " << v[i] << std::endl;
    return os;
}

Vector operator*(double f, const Vector& v) {
    int n = v.size();
    Vector w(n);
    for (int i = 0; i < n; i++)
	w[i] = f*v[i];
    return w;
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

Vector vector_product(const Vector& v1, const Vector& v2) THROW {
    ASSERT(v1.size() == 3 && v2.size() == 3, "v1.size() = " << v1.size() << ", v2.size() = " << v2.size());

    std::vector<double> v(3);
    v[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return Vector(v);
}
