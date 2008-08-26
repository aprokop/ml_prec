#include "prec.h"
#include "include/logger.h"

#include <iomanip>

DEFINE_LOGGER("Prec");

void Prec::solve(Vector& f, Vector& x) THROW {
    Vector _f(f);
    solve(_f, x, 0);
}

void Prec::solve(Vector& f, Vector& x, uint level) THROW {
    Level& li = levels[level];
    uint N = li.N;
    ASSERT(f.size() == N && x.size() == N, "Wrong dimension: N = " << N << ", f = " << f.size() << ", x = " << x.size());

    const std::vector<uint>& dtr   = li.dtr;
    const std::vector<uint>& tr    = li.tr;
    const std::vector<Tail>& tails = li.tails;

    // truncate tails
    for (uint i = 0; i < tails.size(); i++) {
	const Tail& tail = tails[i];
	f[tail[0].index] *= tail[0].a2;

	for (uint j = 1; j < tail.size(); j++) {
	    const TailNode& tn = tail[j];
	    f[tn.index] = tn.a2*f[tn.index] + tn.a3*f[tail[j-1].index];
	}
    }

    if (level < nlevels-1) {
	uint n = levels[level+1].N;

	Vector& f1 = li.f1; 
	for (uint i = 0; i < n; i++) 
	    f1[i] = f[tr[i]];

	Vector& x1 = li.x1; 
	if (ncheb) {
	    // Perform Chebyshev iterations
	    Vector& u0 = li.u0; 
	    Vector& u1 = li.u1;
	    const CSRMatrix& A = levels[level+1].A;

	    double lmin = levels[level+1].lmin;
	    double lmax = levels[level+1].lmax;
	    double eta = (lmax + lmin) / (lmax - lmin);

	    Vector tmp(f1);
	    double alpha, beta;
	    // ===============    STEP 1    ===============
	    alpha = 2/(lmax + lmin);

	    solve(tmp, x1, level+1);
	    x1 *= alpha;

	    // ===============    STEP 2    ===============
	    if (ncheb > 1) {
		u1 = x1;
		alpha = 4/(lmax - lmin) * cheb(eta, 1)/cheb(eta, 2);
		beta  = cheb(eta, 0) / cheb(eta, 2);

		multiply(A, u1, tmp);
		tmp -= f1;
		solve(tmp, x1, level+1);
		for (uint k = 0; k < n; k++) 
		    x1[k] = u1[k] - alpha*x1[k] + beta*u1[k];
	    }

	    // ===============    STEPS 3+    ===============
	    for (uint i = 3; i <= ncheb; i++) {
		// hack to avoid copying and allocating new memory
		u0.swap(x1);
		u1.swap(u0);
		// end hack

		alpha = 4/(lmax - lmin) * cheb(eta, i-1)/cheb(eta, i);
		beta  = cheb(eta, i-2) / cheb(eta, i);

		// x1 = u1 - alpha*solve(A*u1 - f1, level+1) + beta*(u1 - u0);
		multiply(A, u1, tmp);
		tmp -= f1;
		solve(tmp, x1, level+1);
		for (uint k = 0; k < n; k++) 
		    x1[k] = u1[k] - alpha*x1[k] + beta*(u1[k] - u0[k]);
	    }
	} else {
	    solve(f1, x1, level+1);
	}

	for (uint i = 0; i < n; i++)
	    x[tr[i]] = x1[i];

    } else {
	// for last level assert for now that we have only diagonal matrix
	ASSERT(tr.size() == 0, "We must have an \"almost\" diagonal on the coarsest level");
    }

    // solve diagonal part
    for (uint i = 0; i < dtr.size(); i++) {
	uint ii = dtr[i];
	x[ii] = f[ii] / li.aux[ii];
    }

    // restore tails
    for (int i = tails.size()-1; i >= 0; i--) {
	const Tail& tail = tails[i];
	uint tsize = tail.size();
	if (tail.end_type == 'f')
	    x[tail[tsize-1].index] = f[tail[tsize-1].index];
	for (int j = tsize - 2; j >= 0; j--) {
	    const TailNode& tn = tail[j];
	    x[tn.index] = tn.a1*x[tail[j+1].index] + f[tn.index];
	}
    }
}
