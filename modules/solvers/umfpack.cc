#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

#include "umfpack.h"

DEFINE_LOGGER("DirectSolver");

void convert(const CSRMatrix& A, int* ia, int* ja, double* a);
void DirectSolver(const CSRMatrix& A, const Vector& b, Vector& x, void *Symbolic, void *Numeric,
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

    /* FIXME: which parameter should be here? */
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
    umfpack_di_free_numeric(&Numeric);
    umfpack_di_free_symbolic(&Symbolic);

    delete[] ia;
    delete[] ja;
    delete[] a;
}

/* Convert CSR matrix to CSC */
void convert(const CSRMatrix& A, int* ia, int* ja, double* a) {
    const uint n = A.size();

    const uvector<uint>&   Aia = A.get_ia();
    const uvector<uint>&   Aja = A.get_ja();
    const uvector<double>&  Aa = A.get_a();

    std::vector<int> marker(n, 0);

    /* Calculate number of elements in each column */
    for (uint i = 0; i < n; i++)
	for (uint j = Aia[i]; j < Aia[i+1]; j++)
	    marker[Aja[j]]++;

    /*
     * Fill array ia for transposed matrix using marker
     * After cycle marker = B.ia (but shorter by 1 element)
     */
    ia[0] = 0;
    for (uint i = 0; i < n; i++) {
	ia[i+1] = ia[i] + marker[i];
	marker[i] = ia[i];
    }

    /* Fill in ja and a */
    int col, relpos;
    for (uint i = 0; i < n; i++)
	for (uint j = Aia[i]; j < Aia[i+1]; j++) {
	    col = Aja[j];
	    relpos = marker[col];

	    ja[relpos] = i;
	    a[relpos] = Aa[j];

	    marker[col]++;
	}
}
