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

struct SolverStats {
    double t_resid;	    /* Residual calculation time */
    double t_resid_total;
    double t_prec;	    /* Preconditioner inversion time */
    double t_prec_total;
    double t_const;	    /* Preconditioner possible construction time */
    double t_sol;	    /* Time to solve the system */
    uint   niter;	    /* Number of iterations */

    SolverStats() {
	t_resid = t_resid_total = t_prec = t_prec_total = t_const = t_sol = 0;
	niter = 0;
    }
};
std::ostream& operator<<(std::ostream& os, const SolverStats& stats);

void PCGSolver(const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x,
	       double eps = 1e-10, NormType norm_type = NORM_L2, bool silent = false) THROW;
void ChebSolver(const CSRMatrix& A, double lmin, double lmax, const Vector& b, const PrecBase& B, Vector& x,
		double eps = 1e-10, NormType norm_type = NORM_L2, bool silent = false) THROW;
void SimpleSolver(const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x, SolverStats& stats,
		  double eps = 1e-10, NormType norm_type = NORM_L2, bool silent = false) THROW;
void DirectSolver(const CSRMatrix& A, const Vector& b, Vector& x, void *&Symbolic, void *&Numeric,
		  SolverStats& stats, bool silent = false) THROW;

void generate_x0(Vector& x);
double calculate_norm(const Vector& r, const CSRMatrix& A, const PrecBase& B, NormType norm_type);
void check_and_replace_eps(double init_norm, double& eps);

#endif // __SOLVERS_H__
