#include "main.h"
#include "config/config.h"
#include "include/time.h"
#include "modules/prec/prec.h"
#include "modules/solvers/solvers.h"

/* Logger header files */
#include "include/logger.h"
#ifndef NO_LOGGER
#include <log4cxx/propertyconfigurator.h>
#endif

DEFINE_LOGGER("Main");

int main (int argc, char * argv[]) {
    /* Set problem parameters */
    Config cfg;
    if (set_params(argc, argv, cfg)) {
	LLL_INFO("Error while setting parameters, exiting...");
	return 1;
    }

#ifndef NO_LOGGER
    {
	/* Reset logger output file */
	std::ostringstream os;
	std::string regexp_dir = cfg.dir;
	size_t j = regexp_dir.find_first_of('/');
	while(j != std::string::npos) {
	    regexp_dir.replace(j, 1, "\\/");
	    j = regexp_dir.find_first_of('/', j + 2);

	}
	os << "sed 's/log4j.appender.R.File=trace.log/log4j.appender.R.File=" << regexp_dir
		<< "\\/trace.log/' default.properties > log4cxx.properties";
	if (system(os.str().c_str()) != 0)
	    THROW_EXCEPTION("Cannot use sed");
    }

    /* Initialize logger */
    log4cxx::PropertyConfigurator::configure("./log4cxx.properties");
#endif

    LLL_DEBUG("Config parameters:\n" << cfg);

    SPEMesh mesh(cfg.nx, cfg.ny, cfg.nz);

    SkylineMatrix A;
    construct_matrix(cfg, mesh, A);
    // scale_c(A, 0.001);

    Vector b(A.size(), 0.);
    construct_vector(cfg, b);

    if (cfg.dump_data) {
	dump_data(A, b);
	return 0;
    }

#if 0
    // LOG_VAR(A);
    transform(A, b);
    // LOG_VAR(A);
    // exit(0);
#endif

    if (cfg.analysis != ANAL_NONE) {
	analyze(A, b, cfg, cfg.analysis);
	return 0;
    }

    PrecBase * B_ = NULL;
    GlobalStats gstats;

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
	    case AMG_PREC         : B_ = new AMGPrec(A);		break;
	    case DIAG_PREC        : B_ = new DiagPrec(A);		break;
	    case GS_PREC          : B_ = new GSPrec(A);			break;
#ifdef HAVE_UMFPACK
	    case BGS_PREC         : B_ = new BGSPrec(A, cfg.nx*cfg.ny);	break;
	    case RBGS_PREC        : B_ = new RBGSPrec(A, cfg.nx*cfg.ny);break;
#endif
	    case MULTI_SPLIT_PREC : B_ = new MultiSplitPrec(A, cfg);	break;
	    case SYM_SPLIT_PREC   : B_ = new SymPrec(A, cfg);		break;
	    case UH_CHEB_PREC     : B_ = new Prec(A, cfg);		break;
	}
	cfinish = pclock();
	ctimes.push_back(cfinish - cstart);
	LLL_INFO("Construction time : " << ctimes.back());
    }
    PrecBase& B = *B_;
    gstats.t_const = avg_time(ctimes);


    /* Log preconditioner stats */
    if (cfg.prec == UH_CHEB_PREC) {
	Prec& Bcheb = dynamic_cast<Prec&>(B);
	LOG_INFO(Bcheb);
#if 0
	Bcheb.graph_planes("grids.ps", 1, 'z', mesh);
	return 0;
#endif
    } else if (cfg.prec == MULTI_SPLIT_PREC) {
	MultiSplitPrec& Bms = dynamic_cast<MultiSplitPrec&>(B);
	LOG_INFO(Bms);
#if 0
	Bms.graph_planes("grids.ps", 4, 'z', mesh);
	return 0;
#endif
    }

    double eps = 1e-6;

    Vector x(A.size());
    /* =====  Solution phase (preconditioner)  ===== */
    for (uint i = 0; i < cfg.ntests; i++) {
	SolverStats stats;
	switch (cfg.solver) {
	    case PCG_SOLVER : {
		if (cfg.unsym_matrix)
		    LOG_WARN("Applying PCGSolver for unsymmetric matrix");
		PCGSolver(A, b, B, x, eps);
		break;
	    }
	    case CHEB_SOLVER : {
		/*
		 * ChebSolver requires lower and upper bounds on eigenvalues. Thus
		 * it is applicable only to Chebyshev preconditioner
		 */
		Prec& Bcheb = dynamic_cast<Prec&>(B);
		ChebSolver(A, Bcheb.lmin(), Bcheb.lmax(), b, Bcheb, x, eps);
		break;
	    }
	    case SIMPLE_SOLVER : {
		SimpleSolver(A, b, B, x, stats, eps);
		break;
	    }
#ifdef HAVE_UMFPACK
	    case DIRECT_SOLVER : {
		void *Symbolic = NULL, *Numeric = NULL;
		DirectSolver(A, b, x, Symbolic, Numeric, stats);
		break;
	    }
#endif
	}
	stimes.push_back(stats.t_sol);

	gstats.t_prec = stats.t_prec;
	gstats.niter  = stats.niter;
	if (stats.t_const > gstats.t_const)
	    gstats.t_const = stats.t_const;

	LLL_DEBUG(stats);

	LLL_INFO("Solution time : " << stimes.back());
    }

    gstats.t_sol = avg_time(stimes);
    gstats.t_total = gstats.t_const + gstats.t_sol;
    LLL_INFO("Avg construction time: " << gstats.t_const);
    LLL_INFO("Avg solution time    : " << gstats.t_sol);
    LLL_INFO("Avg total time       : " << gstats.t_total);

    delete B_;

    gstats.dump(cfg.dir);

    return 0;
}
