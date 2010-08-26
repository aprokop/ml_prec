#include "main.h"
#include "config/config.h"
#include "include/time.h"
#include "modules/prec/amg/amg_prec.h"
#include "modules/prec/diag/diag.h"
#include "modules/prec/cheb/cheb_prec.h"
#include "modules/prec/relax/rel_prec.h"
#include "modules/prec/sym/sym.h"
#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/solvers/solvers.h"

/* Logger header files */
#include "include/logger.h"
#ifndef NO_LOGGER
#include <log4cxx/propertyconfigurator.h>
#endif

DEFINE_LOGGER("Main");

int main (int argc, char * argv[]) {
#ifndef NO_LOGGER
    /* Initialize logger */
    log4cxx::PropertyConfigurator::configure("./log4cxx.properties");
#endif

    /* Set problem parameters */
    Config cfg;
    if (set_params(argc, argv, cfg)) {
	LLL_INFO("Error while setting parameters, exiting...");
	return 1;
    }
    LLL_DEBUG("Config parameters:\n" << cfg);

    SPEMesh mesh(cfg.nx, cfg.ny, cfg.nz);

    SkylineMatrix A;
    construct_matrix(cfg, mesh, A);

    Vector b(A.size(), 0.);
    construct_vector(cfg, b);

    if (cfg.dump_data) {
	dump_data(A, b);
	return 0;
    }

    PrecBase * B_ = NULL;

    std::vector<double> ctimes, stimes;

    /* =====  Construction phase (preconditioner)  ===== */
    /* Timers */
    double cstart, cfinish, sstart, sfinish;
    for (uint i = 0; i < cfg.ntests; i++) {
	if (B_) {
	    delete B_;
	    B_ = NULL;
	}
	cstart = cfinish = 0;
	cstart = pclock();
	switch (cfg.prec) {
	    case UH_CHEB_PREC     : B_ = new Prec(A, cfg);		break;
	    case AMG_PREC         : B_ = new AMGPrec(A);		break;
	    case DIAG_PREC        : B_ = new DiagPrec(A);		break;
	    case SYM_SPLIT_PREC   : B_ = new SymPrec(A, cfg);		break;
	    case MULTI_SPLIT_PREC : B_ = new MultiSplitPrec(A, cfg);	break;
	}
	cfinish = pclock();
	ctimes.push_back(cfinish - cstart);
	LLL_INFO("Construction time : " << ctimes.back());
    }
    PrecBase& B = *B_;

    /* Log preconditioner stats */
    if (cfg.prec == UH_CHEB_PREC) {
	Prec& Bcheb = static_cast<Prec&>(B);
	LOG_INFO(Bcheb);
#if 0
	Bcheb.graph_planes("grids.ps", 1, 'z', mesh);
	return 0;
#endif
    } else if (cfg.prec == MULTI_SPLIT_PREC) {
	MultiSplitPrec& Bms = static_cast<MultiSplitPrec&>(B);
	LOG_INFO(Bms);
#if 0
	Bms.graph_planes("grids.ps", 1, 'z', mesh);
	return 0;
#endif
    }

    double eps = 1e-6;

    Vector x(A.size());
    /* =====  Solution phase (preconditioner)  ===== */
    for (uint i = 0; i < cfg.ntests; i++) {
	sstart = pclock();
	if (!cfg.unsym_matrix) {
	    /* The matrix is symmetric */
	    if (cfg.solver == PCG_SOLVER)
		PCGSolver(A, b, B, x, eps);
	    else {
		/* ChebSolver requires lower and upper bounds on eigenvalues. Thus
		 * it is applicable only to Chebyshev preconditioner */
		Prec& Bcheb = static_cast<Prec&>(B);
		ChebSolver(A, Bcheb.lmin(), Bcheb.lmax(), b, Bcheb, x, eps);
	    }
	} else {
	    /*
	     * The matrix is unsymmetric
	     * Works only with a simple solver (exterior) for now. In the future,
	     * we could use them with GMRES of BiCGStab
	     */
	    SimpleSolver(A, b, B, x, eps);
	}
	sfinish = pclock();
	stimes.push_back(sfinish - sstart);

	LLL_INFO("Solution time : " << stimes.back());
    }

    double ctime = avg_time(ctimes);
    double stime = avg_time(stimes);
    LLL_INFO("Avg construction time: " << ctime);
    LLL_INFO("Avg solution time    : " << stime);
    LLL_INFO("Avg total time       : " << ctime + stime);

    delete B_;

    return 0;
}
