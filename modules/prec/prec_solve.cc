#include "prec.h"
#include "include/logger.h"

#include <iomanip>

DEFINE_LOGGER("Prec");

clock_t diag_time = 0;
clock_t mult_time = 0;
clock_t restr_f   = 0;
clock_t prol_x    = 0;
void Prec::solve(Vector& f, Vector& x) THROW {
    diag_time = 0;
    mult_time = 0;
    restr_f   = 0;
    prol_x    = 0;

    solve(f, x, 0);

#if 1
    LOG_DEBUG("Diagonal time = " << std::setprecision(3) << double(diag_time)/CLOCKS_PER_SEC);
    LOG_DEBUG("Multipli time = " << std::setprecision(3) << double(mult_time)/CLOCKS_PER_SEC);
    LOG_DEBUG("Restrict time = " << std::setprecision(3) << double(restr_f)/CLOCKS_PER_SEC);
    LOG_DEBUG("Prolong  time = " << std::setprecision(3) << double(prol_x)/CLOCKS_PER_SEC);
#endif
}

void Prec::solve(Vector f, Vector& x, uint level) THROW {
    Level& li = levels[level];
    uint N = li.N;
    ASSERT(f.size() == N && x.size() == N, "Wrong dimension: N = " << N << ", f = " << f.size() << ", x = " << x.size());

    clock_t delta, gdelta = clock();

    const std::vector<uint>& dtr = li.dtr;
    const std::vector<uint>& tr  = li.tr;
    std::vector<Tail>& tails     = li.tails;

    for (uint i = 0; i < tails.size(); i++) {
	Tail& tail = tails[i];
	f[tail[0].index] *= tail[0].a2;

	for (uint j = 1; j < tail.size(); j++) {
	    TailNode& tn = tail[j];
	    f[tn.index] = tn.a2*f[tn.index] + tn.a3*f[tail[j-1].index];
	}
    }

    if (level < nlevels-1) {
	uint n = levels[level+1].N;

	Vector& f1 = li.f1; 
	delta = clock();
	for (uint i = 0; i < n; i++) 
	    f1[i] = f[tr[i]];
	restr_f += clock() - delta;

	Vector& x1 = li.x1; 
	if (ncheb) {
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
	    if (ncheb > 1) {
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
	    for (uint i = 3; i <= ncheb; i++) {
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
	ASSERT(tr.size() == 0, "We must have an \"almost\" diagonal on the coarsest level");
    }

    // solve diagonal part
    delta = clock();
    for (uint i = 0; i < dtr.size(); i++) {
	uint ii = dtr[i];
	x[ii] = f[ii] / li.aux[ii];
    }
    diag_time += clock() - delta;

    // restore tails
    for (int i = tails.size()-1; i >= 0; i--) {
	const Tail& tail = tails[i];
	uint tsize = tail.size();
	if (tail.end_is_local == false) 
	    x[tail[tsize-1].index] = f[tail[tsize-1].index];
	for (int j = tsize - 2; j >= 0; j--) {
	    const TailNode& tn = tail[j];
	    x[tn.index] = tn.a1*x[tail[j+1].index] + f[tn.index];
	}
    }
}
