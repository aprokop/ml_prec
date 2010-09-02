#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"
#include "modules/prec/prec_base.h"

enum NormType {
    NORM_L2,
    NORM_A,
    NORM_B_1
};

void PCGSolver(const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x,
	       double eps = 1e-10, NormType norm_type = NORM_L2, bool silent = false) THROW;
void ChebSolver(const CSRMatrix& A, double lmin, double lmax, const Vector& b, const PrecBase& B, Vector& x,
		double eps = 1e-10, NormType norm_type = NORM_L2, bool silent = false) THROW;
void SimpleSolver(const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x,
		  double eps = 1e-10, NormType norm_type = NORM_L2, bool silent = false) THROW;

void generate_x0(Vector& x);
double calculate_norm(const Vector& r, const CSRMatrix& A, const PrecBase& B, NormType norm_type);
void check_and_replace_eps(double init_norm, double& eps);

#endif // __SOLVERS_H__
