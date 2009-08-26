#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "config/config.h"
#include "include/define.h"
#include "include/exception.h"
#include "include/uvector.h"

typedef uvector<double> Vector;

bool is_nan(const Vector& v);

void dump(const std::string& filename, bool ascii = false);
void load(const std::string& filename, bool ascii = false);

void   daxpy(double alpha, const Vector& x, Vector& y);
void   dscal(double alpha, Vector& x);
double ddot(const Vector& x, const Vector& y);
double dnrm2(const Vector& x);

/* Simple variants */
void   daxpy(double alpha, const double* x, double* y, uint n);
void   dscal(double alpha, double* x, uint n);
double ddot(const double* x, const double* y, uint n);

Vector vector_product(const Vector& v1, const Vector& v2) THROW;
std::ostream& operator<<(std::ostream& os, const Vector& v);

#endif // #ifndef __VECTOR_H__
