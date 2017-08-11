#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "config/config.h"
#include "include/define.h"
#include "include/exception.h"
#include "include/uvector.h"

typedef uvector<double> Vector;

bool is_nan(const Vector& v);

void   daxpy(double alpha, const Vector& x, Vector& y);
void   dscal(double alpha, Vector& x);
double ddot(const Vector& x, const Vector& y);
double dnrm2(const Vector& x);

/* Simple variants */
void   daxpy(double alpha, const double* x, double* y, uint n);
void   dscal(double alpha, double* x, uint n);
double ddot(const double* x, const double* y, uint n);
double dnrm2(const double* x, uint n);

Vector vector_product(const Vector& v1, const Vector& v2);

enum DumpType {
    ASCII,
    BINARY,
    HYPRE,
    MATRIX_MARKET
};

std::ostream& operator<<(std::ostream& os, const Vector& v);
void dump(const std::string& filename, const Vector& v, DumpType type);
void load(Vector& v, const std::string& filename, DumpType type);

#endif // #ifndef __VECTOR_H__
