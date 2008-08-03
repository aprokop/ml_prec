#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "config/config.h"
#include "include/exception.h"

#include <vector>
#include <tr1/memory>

class Vector {
private:
    typedef std::vector<double> data_type;
    std::tr1::shared_ptr<data_type> data;
    uint n;

    void check_index(uint i) const THROW {
	ASSERT(data->size() == n, "data->size() = " << data->size() << ", n = " << n);
	ASSERT(i < n, "Index is out of bundaries: i = " << i << ", n = " << n);
    }
    void set_unique(void) {
	if (!data.unique())
	    data.reset(new data_type(*data));
    }

public:
    Vector(uint _n = 0);
    Vector(const std::vector<double>& v);
    Vector(const Vector& v);
    ~Vector() {}

    uint size() const;

    // std::vector functions
    void resize(uint new_size, double v = 0.);
    void reserve(uint size);

    const double& operator[](uint index) const THROW {
	check_index(index);
	return (*data)[index];
    }
    double& operator[](uint index) THROW;

    bool operator==(const Vector& v) const THROW;
    bool operator!=(const Vector& v) const THROW;

    const Vector& operator+() const;
    Vector operator-() const;
    Vector operator+(const Vector& v) const THROW;
    Vector operator-(const Vector& v) const THROW;
    Vector operator*(double f) const;
    Vector operator/(double f) const THROW;

    const Vector& operator=(const Vector& v);
    void copy(const Vector& v);

    const Vector& operator+=(const Vector& v) THROW;
    const Vector& operator-=(const Vector& v) THROW;
    const Vector& operator*=(double f);
    const Vector& operator/=(double f) THROW;

    double norm_1() const;
    double norm_2() const;
    double norm_inf() const;
};

std::ostream& operator<<(std::ostream& os, const Vector& v);
Vector operator*(double f, const Vector& v);

double scalar_product(const Vector& v1, const Vector& v2) THROW;
Vector vector_product(const Vector& v1, const Vector& v2) THROW;

#endif // #ifndef __VECTOR_H__
