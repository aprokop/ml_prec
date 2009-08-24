#include "include/logger.h"
#include "cheb_prec.h"
#include "modules/common/common.h"

#include <iomanip>

DEFINE_LOGGER("Prec");

/* NOTE: we'd like to pass f as const but due to having tails it is simplier to do otherwise */
static Vector _f;
void Prec::solve(Vector& f, Vector& x) THROW {
    _f = f;
    solve(0, _f, x);
}

void Prec::solve(uint level, Vector& f, Vector& x) const THROW {
    const Level& li = levels[level];
    const Level& ln = levels[level+1];

    uint N = li.N;

    const uvector<uint>& dtr = li.dtr;
    const uvector<uint>& tr  = li.tr;

    const std::vector<Tail>& tails  = li.tails;
    if (use_tails) {
	/* Truncate tails */
	for (uint i = 0; i < tails.size(); i++) {
	    const Tail& tail = tails[i];
	    f[tail[0].index] *= tail[0].a2;

	    for (uint j = 1; j < tail.size(); j++) {
		const TailNode& tn = tail[j];
		f[tn.index] = tn.a2*f[tn.index] + tn.a3*f[tail[j-1].index];
	    }
	}
    }

    if (level < nlevels-1) {
	uint n = ln.N;
	uint ncheb = li.ncheb;

	Vector& f1 = li.f1; 
	for (uint i = 0; i < n; i++) 
	    f1[i] = f[tr[i]];


	// Perform Chebyshev iterations
	const CSRMatrix& A = levels[level+1].A;

	double lmin = levels[level+1].lmin;
	double lmax = levels[level+1].lmax;
	double eta  = (lmax + lmin) / (lmax - lmin);

	Vector& tmp = li.tmp;
	tmp = f1;
	double alpha, beta;

	Vector& x1 = li.x1;
	Vector& u1 = li.u1;
	// ===============    STEP 1    ===============
	alpha = 2/(lmax + lmin);

	if (ncheb > 1) {
	    solve(level+1, tmp, u1);
	    dscal(alpha, u1);
	} else {
	    solve(level+1, tmp, x1);
	    dscal(alpha, x1);
	}

	// ===============    STEP 2    ===============
	if (ncheb > 1) {
	    /* x1 = (1 + beta)*u1 - alpha*solve(A*u1 - f1, level+1) */
	    alpha = 4/(lmax - lmin) * cheb(eta, 1)/cheb(eta, 2);
	    beta  = cheb(eta, 0) / cheb(eta, 2);


	    multiply(A, u1, tmp, 's');
	    daxpy(-1, f1, tmp);
	    solve(level+1, tmp, x1);
	    dscal(-alpha, x1);
	    daxpy(1 + beta, u1, x1);
	}

	// ===============    STEPS 3+    ===============
	Vector& u0 = li.u0;
	for (uint i = 3; i <= ncheb; i++) {
	    /* Use little hack to avoid copying and allocating new memory */
	    u0.swap(x1);
	    u1.swap(u0);
	    /* end hack */

	    alpha = 4/(lmax - lmin) * cheb(eta, i-1)/cheb(eta, i);
	    beta  = cheb(eta, i-2) / cheb(eta, i);

	    /* x1 = u1 - alpha*solve(A*u1 - f1, level+1) + beta*(u1 - u0) */
	    multiply(A, u1, tmp, 's');
	    daxpy(-1, f1, tmp);
	    solve(level+1, tmp, x1);
	    for (uint k = 0; k < n; k++) 
		x1[k] = u1[k] - alpha*x1[k] + beta*(u1[k] - u0[k]);
	}

	/* Copy the result of Chebyshev iterations to x */
	for (uint i = 0; i < n; i++)
	    x[tr[i]] = x1[i];
    } 

    /* Solve diagonal subsystem */
    for (uint _i = 0; _i < dtr.size(); _i++) {
	uint i = dtr[_i];
	x[i] = f[i] / li.aux[i];
    }

    if (use_tails) {
	/* Restore tails in *backwards* order */
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
}
