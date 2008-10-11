#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"
#include "modules/prec/prec_base.h"

Vector PCGSolver(const CSRMatrix& A, const Vector& b, PrecBase& B, double eps = 1e-10) THROW;
Vector ChebSolver(const CSRMatrix& A, double lmin, double lmax, const Vector& b, PrecBase& B, double eps) THROW;

#endif // __SOLVERS_H__
