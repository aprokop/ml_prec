#include "solvers.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

DEFINE_LOGGER("PCG");

void PCGSolver(const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x, SolverStats& stats,
	       double eps, NormType norm_type, bool silent) THROW {
    double  gtime = pclock();
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    uint   n = b.size();
    Vector r(n), z(n);
    Vector p(n), Ap(n);
    double alpha, beta;
    double app, rz0, rz1;
    double norm, init_norm;

    double  mult = 0,  inv = 0,  cstr = 0, delta;
    int	   nmult = 0, ninv = 0;

    generate_x0(x);
    residual(A, b, x, r);
    norm = init_norm = calculate_norm(r, A, B, norm_type);

    delta = pclock();
    B.solve(r, z);
    cstr = pclock() - delta;

    p = z;
    rz0 = ddot(r, z);

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
	multiply(A, p, Ap);
	mult += pclock() - delta;
	nmult++;

	app = ddot(Ap, p);
	if (app <= 0) {
	    LOG_ERROR("<Ap, p> <= 0");
	    THROW_EXCEPTION("<Ap, p> = " << app << ": |Ap| = " << dnrm2(Ap) << ", |p| = " << dnrm2(p));
	}
	alpha = rz0 / app;

	daxpy(alpha, p, x);
	daxpy(-alpha, Ap, r);

	norm = calculate_norm(r, A, B, norm_type);

	delta = pclock();
	B.solve(r, z);
	inv += pclock() - delta;
	ninv++;

	rz1 = ddot(r, z);
	beta = rz1 / rz0;

	dscal(beta, p);
	daxpy(1, z, p);

	rz0 = rz1;
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
