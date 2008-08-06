#include "prec.h"
#include "include/logger.h"

#include <iomanip>

DEFINE_LOGGER("Prec");

clock_t diag_time = 0;
clock_t mult_time = 0;
clock_t restr_f   = 0;
clock_t prol_x    = 0;
clock_t c_time	  = 0;
void Prec::solve(Vector& f, Vector& x) THROW {
    diag_time = 0;
    mult_time = 0;
    restr_f   = 0;
    prol_x    = 0;
    c_time    = 0;

    solve(f, x, 0);

#if 0
    LOG_DEBUG("Diagonal time = " << std::setprecision(3) << double(diag_time)/CLOCKS_PER_SEC);
    LOG_DEBUG("Multipli time = " << std::setprecision(3) << double(mult_time)/CLOCKS_PER_SEC);
    LOG_DEBUG("Restrict time = " << std::setprecision(3) << double(restr_f)/CLOCKS_PER_SEC);
    LOG_DEBUG("Prolong  time = " << std::setprecision(3) << double(prol_x)/CLOCKS_PER_SEC);
    LOG_DEBUG("Coarse   time = " << std::setprecision(3) << double(c_time)/CLOCKS_PER_SEC);
#endif
}

void Prec::solve(const Vector& f, Vector& x, uint level) THROW {
    Level& li = levels[level];
    uint N = li.N;
    ASSERT(f.size() == N && x.size() == N, "Wrong dimension: N = " << N << ", f = " << f.size() << ", x = " << x.size());

    clock_t delta, gdelta = clock();
    const std::vector<uint>& dtr = li.dtr;
    if (level < nlevels-1) {
	const std::vector<uint>& tr = li.tr;

	uint n = levels[level+1].N;

	Vector& f1 = li.f1; 
	delta = clock();
	for (uint i = 0; i < n; i++) 
	    f1[i] = f[tr[i]];
	restr_f += clock() - delta;

	Vector& x1 = li.x1; 
	if (li.ncheb) {
	    // Perform Chebyshev iterations
	    Vector& u0 = li.u0; 
	    Vector& u1 = li.u1;
	    const CSRMatrix& A = levels[level+1].A;

	    double lmin = levels[level+1].lmin;
	    double lmax = levels[level+1].lmax;
	    double eta = (lmax + lmin) / (lmax - lmin);

	    Vector tmp(n);
	    double alpha, beta;
	    // ===============    STEP 1    ===============
	    alpha = 2/(lmax + lmin);

	    solve(f1, x1, level+1);
	    x1 *= alpha;

	    // ===============    STEP 2    ===============
	    if (li.ncheb > 1) {
		u1.copy(x1);
		alpha = 4/(lmax - lmin) * cheb(eta, 1)/cheb(eta, 2);
		beta  = cheb(eta, 0) / cheb(eta, 2);

		delta = clock();
		multiply(A, u1, tmp);
		mult_time += clock() - delta;
		tmp -= f1;
		solve(tmp, x1, level+1);
		for (uint k = 0; k < n; k++) 
		    x1[k] = u1[k] - alpha*x1[k] + beta*u1[k];
	    }

	    // ===============    STEPS 3+    ===============
	    for (uint i = 3; i <= li.ncheb; i++) {
		// hack to avoid copying and allocating new memory
		u0.swap(x1);
		u1.swap(u0);
		// end hack

		alpha = 4/(lmax - lmin) * cheb(eta, i-1)/cheb(eta, i);
		beta  = cheb(eta, i-2) / cheb(eta, i);

		// x1 = u1 - alpha*solve(A*u1 - f1, level+1) + beta*(u1 - u0);
		delta = clock();
		multiply(A, u1, tmp);
		mult_time += clock() - delta;
		tmp -= f1;
		solve(tmp, x1, level+1);
		for (uint k = 0; k < n; k++) 
		    x1[k] = u1[k] - alpha*x1[k] + beta*(u1[k] - u0[k]);
	    }
	} else {
	    solve(f1, x1, level+1);
	}

	delta = clock();
	for (uint i = 0; i < n; i++)
	    x[tr[i]] = x1[i];
	prol_x += clock() - delta;

    } else {
	// for last level assert for now that we have only diagonal matrix
	ASSERT(dtr.size() == N, "We must have a diagonal on the coarsest level");
    }

    // solve diagonal part
    delta = clock();
    for (uint i = 0; i < dtr.size(); i++) 
	x[dtr[i]] = f[dtr[i]] / c;
    diag_time += clock() - delta;

    if (level > 6)
	c_time += clock() - gdelta;
}

