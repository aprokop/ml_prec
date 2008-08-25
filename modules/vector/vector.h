#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "config/config.h"
#include "include/exception.h"

#include <vector>
#include <tr1/memory>

class Vector {
private:
    std::vector<double> data;

public:
    Vector(uint n = 0) : data(n) { }
    ~Vector() { }
    const Vector& operator=(const Vector& v);

    const double& operator[](uint i) const THROW {
	ASSERT(i < data.size(), "Index is out of bundaries: i = " << i << ", n = " << data.size());
	return data[i];
    }
    double& operator[](uint i) THROW {
	ASSERT(i < data.size(), "Index is out of bundaries: i = " << i << ", n = " << data.size());
	return data[i];
    }

    uint size() const {
	return data.size();
    }
    void resize(uint n) {
	data.resize(n, 0.);
    }
    void swap(Vector& v) {
	data.swap(v.data);
    }

    const Vector& operator+=(const Vector& v) THROW;
    const Vector& operator-=(const Vector& v) THROW;
    const Vector& operator*=(double f);
    const Vector& operator/=(double f) THROW;

    double norm_2() const;
};

double scalar_product(const Vector& v1, const Vector& v2) THROW;

#endif // #ifndef __VECTOR_H__
