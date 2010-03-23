#include "include/logger.h"
#include "multi_split_prec.h"
#include "modules/common/common.h"

#include <iomanip>

DEFINE_LOGGER("MultiSplitPrec");

void MultiSplitPrec::truncate_tails(uint level, Vector& f) const THROW {
    const std::vector<Tail>& tails = levels[level].tails;

    for (uint i = 0; i < tails.size(); i++) {
	const Tail& tail = tails[i];
	f[tail[0].index] *= tail[0].a2;

	for (uint j = 1; j < tail.size(); j++) {
	    const TailNode& tn = tail[j];
	    f[tn.index] = tn.a2*f[tn.index] + tn.a3*f[tail[j-1].index];
	}
    }
}

void MultiSplitPrec::restore_tails(uint level, const Vector& f, Vector& x) const THROW {
    const std::vector<Tail>& tails = levels[level].tails;

    for (int i = tails.size()-1; i >= 0; i--) {
	const Tail& tail = tails[i];
	if (tail.end_type == 'f')
	    x[tail.back().index] = f[tail.back().index];
	for (int j = tail.size() - 2; j >= 0; j--) {
	    const TailNode& tn = tail[j];
	    x[tn.index] = tn.a1*x[tail[j+1].index] + f[tn.index];
	}
    }
}

void MultiSplitPrec::solve_diagonal(uint level, const Vector& f, Vector& x) const THROW {
    const Level& li                = levels[level];
    const uvector<uint>& dtr       = li.dtr;
    const uvector<double>& dtr_val = li.dtr_val;

    for (uint i_ = 0; i_ < dtr.size(); i_++) {
	uint i = dtr[i_];
	x[i] = f[i] * dtr_val[i_];
    }
}

void MultiSplitPrec::solve(Vector& f, Vector& x) THROW {
    solve(0, f, x);
}

void MultiSplitPrec::solve(uint level, const Vector& f, Vector& x) const THROW {
    const Level& li = levels[level];

    uint N = li.N;
    uint n = (level < nlevels-1) ? levels[level+1].N : 0;

    const uvector<uint>& tr  = li.tr;
    uint niter = li.niter;

    const CSRMatrix& A = levels[level].A;

    Vector& r  = li.r;
    Vector& u0 = li.u0;
    Vector& u1 = li.u1;
    Vector& r1 = li.r1;

    /* Copy f into modifiable vector */
    r = f;

    /* ===============    STEP 1    =============== */
    /*
     * x1 = x0 + solve(f1 - A*x0)
     * But x0 = 0, so x1 = solve(f1)
     */
    truncate_tails(level, r);
    if (level < nlevels-1) {
	for (uint i = 0; i < n; i++)
	    r1[i] = r[tr[i]];

	solve(level+1, r1, u1);

	for (uint i = 0; i < n; i++)
	    x[tr[i]] = u1[i];
    }
    solve_diagonal(level, r, x);
    restore_tails(level, r, x);

    /* ===============    STEP 2+    =============== */
    /* x^{k+1} = x^k + solve(f1 - A*x^k) */
    for (uint i = 2; i < niter; i++) {
	residual(A, f, x, r);

	truncate_tails(level, r);
	if (level < nlevels-1) {
	    for (uint i = 0; i < n; i++)
		r1[i] = r[tr[i]];

	    solve(level+1, r1, u1);

	    for (uint i = 0; i < n; i++)
		u0[tr[i]] = u1[i];
	}

	solve_diagonal(level, r, u0);
	restore_tails(level, r, u0);

	daxpy(1., u0, x);
    }
}
