#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

#include "umfpack.h"

DEFINE_LOGGER("DirectSolver");

void DirectSolver(const CSRMatrix& A, const Vector& b, Vector& x, void * &Symbolic, void * &Numeric,
                  SolverStats& stats, bool silent) THROW {
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    const int n = A.size();
    const int nnz = A.nnz();

    /* Convert matrix to the CSC format */
    int *ia   = new int[n+1];
    int *ja   = new int[nnz];
    double *a = new double[nnz];
    convert(A, ia, ja, a);

    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];

    /* Set the default control parameters in the Control array */
    umfpack_di_defaults(Control);

    Control[UMFPACK_PRL] = 4;

    if (Symbolic == NULL || Numeric == NULL) {
        /* Perform symbolic factorization */
        umfpack_di_symbolic(n, n, ia, ja, a, &Symbolic, 0, Info);

        /* Perform numeric factorization */
        umfpack_di_numeric(ia, ja, a, Symbolic, &Numeric, 0, Info);
    }

    /* Solve the system */
    umfpack_di_solve(UMFPACK_A, ia, ja, a, &x[0], &b[0], Numeric, 0, Info);

    /* Free memory */
    // umfpack_di_free_numeric(Numeric);
    // umfpack_di_free_symbolic(Symbolic);

    delete[] ia;
    delete[] ja;
    delete[] a;
}

void DirectSolver_free(void * &Symbolic, void * &Numeric) {
    if (Symbolic) { umfpack_di_free_symbolic(&Symbolic); Symbolic = NULL; }
    if (Numeric)  { umfpack_di_free_numeric(&Numeric);   Numeric = NULL; }
}
