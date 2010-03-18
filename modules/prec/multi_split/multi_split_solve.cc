#include "include/logger.h"
#include "multi_split_prec.h"
#include "modules/common/common.h"

#include <iomanip>

DEFINE_LOGGER("MultiSplitPrec");

/* NOTE: we'd like to pass f as const but due to having tails it is simplier to do otherwise */
static Vector f_;
void MultiSplitPrec::solve(Vector& f, Vector& x) THROW {
    f_ = f;
    solve(0, f_, x);
}

void MultiSplitPrec::solve(uint level, Vector& f, Vector& x) const THROW {
    const Level& li = levels[level];
    const Level& ln = levels[level+1];

    uint N = li.N;

    const uvector<uint>& tr  = li.tr;

    if (level < nlevels-1) {
	uint n = ln.N;
	uint niter = li.niter;

	const CSRMatrix& A = levels[level+1].A;

	Vector& f1 = li.f1;
	Vector& x1 = li.x1;
	Vector& u0 = li.u0;
	Vector& r  = li.r;

	/* Construct rhs */
	for (uint i = 0; i < n; i++)
	    f1[i] = f[tr[i]];

	/* ===============    STEP 1    =============== */
	/*
	 * x1 = x0 + solve(f1 - A*x0)
	 * But x0 = 0, so x1 = solve(f1)
	 */
	solve(level+1, f1, x1);

	/* ===============    STEP 2+    =============== */
	/* x^{k+1} = x^k + solve(f1 - A*x^k) */
	for (uint i = 2; i < niter; i++) {
	    residual(A, f1, x1, r);
	    solve(level+1, r, u0);
	    daxpy(1., u0, x1);
	}

	/* Copy the result of iterations to x */
	for (uint i = 0; i < n; i++)
	    x[tr[i]] = x1[i];
    }

    /* Solve diagonal subsystem */
    const uvector<uint>& dtr = li.dtr;
    const uvector<double>& dtr_val = li.dtr_val;
    for (uint i_ = 0; i_ < dtr.size(); i_++) {
	uint i = dtr[i_];
	x[i] = f[i] * dtr_val[i_];
    }
}
