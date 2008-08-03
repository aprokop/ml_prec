#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include "modules/matrix/matrix.h"
#include "modules/prec/prec.h"

Vector CG(const MatrixInterface& A, const Vector& b, double eps = 1e-10) THROW;
Vector PCG(const MatrixInterface& A, const Vector& b, PrecBase& B, double eps = 1e-10) THROW;

#endif // __SOLVERS_H__
