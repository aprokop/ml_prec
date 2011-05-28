#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("SimpleSolver");

void SimpleSolver(const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x, SolverStats& stats,
		  double eps, NormType norm_type, bool silent) THROW {
    double  gtime = pclock();
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());
    if (norm_type != NORM_L2)
	THROW_EXCEPTION("One must only use L2 norm in simple solver");

    uint   n = b.size();
    Vector r(n), z(n);
    double norm, init_norm;

    double  mult = 0,  inv = 0,  cstr = 0, delta;
    int	   nmult = 0, ninv = 0;

#if 1
    generate_x0(x);
#else
    memset(&x[0], 0, x.size()*sizeof(double));
#endif
    residual(A, b, x, r);

    norm = init_norm = calculate_norm(r, A, B, norm_type);

    int niter = 0;
#ifdef ABSOLUTE_NORM
    init_norm = 1;
#else
    check_and_replace_eps(init_norm, eps);
#endif

    while (norm/init_norm > eps) {
	if (!silent)
	    LOG_DEBUG("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

	delta = pclock();
	B.solve(r, z);
	if (niter) {
	    inv += pclock() - delta;
	    ninv++;
	} else {
	    cstr = pclock() - delta;
	}

	daxpy(1., z, x);

	delta = pclock();
	residual(A, b, x, r);
	mult += pclock() - delta;
	nmult++;

	norm = calculate_norm(r, A, B, norm_type);

	niter++;
    }
    gtime = pclock() - gtime;

    if (ninv) {
	stats.t_resid = mult/nmult;
	stats.t_prec  = inv/ninv;
	stats.t_const = 0;

	stats.niter         = niter;
	stats.t_resid_total = mult;
	stats.t_prec_total  = inv;
	stats.t_sol         = gtime;

	double cstr_pos = cstr - stats.t_prec;
	if (cstr_pos > 1.05*stats.t_prec) {
	    stats.t_const  = cstr_pos;
	    stats.t_sol   -= cstr_pos;
	}
    } else {
	THROW_EXCEPTION("No need to iterate");
    }

    if (!silent)
	LLL_INFO("#" << stats.niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
    else
	LOG_INFO("#" << stats.niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);
}
