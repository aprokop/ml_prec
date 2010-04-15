#include "include/logger.h"
#include "cheb_prec.h"
#include "modules/common/common.h"

#include <iomanip>

DEFINE_LOGGER("Prec");

void Prec::solve(Vector& f, Vector& x) THROW {
    solve(0, f, x);
}

void Prec::solve(uint level, const Vector& f, Vector& x) const THROW {
    const Level& li = levels[level];
    const Level& ln = levels[level+1];

    uint N = li.N;
    uint M = li.M;

    const uvector<uint>& map  = li.map;

    const SkylineMatrix& U = li.U;
    const CSRMatrix&     L = li.L;

    /* Solve L*w = f */
    Vector& w = li.w;
    for (uint i = 0; i < N; i++) { /* i is a permuted index */
	w[i] = f[map[i]];
	for (uint j_ = L.ia[i]; j_ < L.ia[i+1]; j_++)
	    w[i] -= L.a[j_] * w[L.ja[j_]];
    }

    Vector& x2 = li.x2;
    if (level < nlevels-1) {
	uint n = ln.N;
	uint ncheb = li.ncheb;

	uvector<double> F(n);
	memcpy(&F[0], &w[M], n*sizeof(double));

	/* Perform Chebyshev iterations */
	const CSRMatrix& A = levels[level+1].A;

	double lmin = levels[level+1].lmin;
	double lmax = levels[level+1].lmax;
	double eta  = (lmax + lmin) / (lmax - lmin);

	Vector& tmp = li.tmp;
	tmp = F;
	double alpha, beta;

	Vector& u1 = li.u1;
	// ===============    STEP 1    ===============
	alpha = 2/(lmax + lmin);

	if (ncheb > 1) {
	    solve(level+1, tmp, u1);
	    dscal(alpha, u1);
	} else {
	    solve(level+1, tmp, x2);
	    dscal(alpha, x2);
	}

	// ===============    STEP 2    ===============
	if (ncheb > 1) {
	    /* x2 = (1 + beta)*u1 - alpha*solve(A*u1 - F, level+1) */
	    alpha = 4/(lmax - lmin) * cheb(eta, 1)/cheb(eta, 2);
	    beta  = cheb(eta, 0) / cheb(eta, 2);

	    multiply(A, u1, tmp, 's');
	    daxpy(-1, F, tmp);
	    solve(level+1, tmp, x2);
	    dscal(-alpha, x2);
	    daxpy(1 + beta, u1, x2);
	}

	// ===============    STEPS 3+    ===============
	Vector& u0 = li.u0;
	for (uint i = 3; i <= ncheb; i++) {
	    /* Use little hack to avoid copying and allocating new memory */
	    u0.swap(x2);
	    u1.swap(u0);
	    /* end hack */

	    alpha = 4/(lmax - lmin) * cheb(eta, i-1)/cheb(eta, i);
	    beta  = cheb(eta, i-2) / cheb(eta, i);

	    /* x2 = u1 - alpha*solve(A*u1 - F, level+1) + beta*(u1 - u0) */
	    multiply(A, u1, tmp, 's');
	    daxpy(-1, F, tmp);
	    solve(level+1, tmp, x2);
	    for (uint k = 0; k < n; k++)
		x2[k] = u1[k] - alpha*x2[k] + beta*(u1[k] - u0[k]);
	}

	for (uint i = M; i < N; i++)
	    x[map[i]] = x2[i-M];
    }

    for (int i = M-1; i >= 0; i--) {
	double z = w[i];

	uint j_;
	for (j_ = U.ia[i+1]-1; j_ > U.ia[i]; j_--) {
	    z -= U.a[j_] * x[map[U.ja[j_]]];
	}

	x[map[i]] = z / U.a[U.ia[i]];
    }
}
