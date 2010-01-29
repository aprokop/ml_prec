#include "main.h"
#include "config/config.h"
#include "include/time.h"
#include "modules/mesh/mesh.h"
#include "modules/prec/amg/amg_prec.h"
#include "modules/prec/diag/diag.h"
#include "modules/prec/cheb/cheb_prec.h"
#include "modules/prec/relax/rel_prec.h"
#include "modules/solvers/solvers.h"

// logger
#include "include/logger.h"
#ifndef NO_LOGGER
#include <log4cxx/propertyconfigurator.h>
#endif

DEFINE_LOGGER("Main");

int main (int argc, char * argv[]) {
#ifndef NO_LOGGER
    // Initialize logger
    log4cxx::PropertyConfigurator::configure("./log4cxx.properties");
#endif

    Config cfg;
    if (set_params(argc, argv, cfg)) {
	LOG_INFO("Error while setting parameters, exiting...");
	std::cout << "Error while setting parameters, exiting..." << std::endl;
	return 1;
    }

    SkylineMatrix A;

    if (cfg.matrix.empty()) {
	SPEMesh mesh;
	mesh.construct_matrix(A, cfg.c);
    } else {
	/* Whether we read matrix in CSR format (transform = true) or already in Skyline (false) */
	bool transform = true;
	A.load(cfg.matrix, transform);
    }

    std::cout << cfg << std::endl;
    LOG_DEBUG("Config parameters: " << cfg);

    PrecBase * B_ = NULL;

    std::vector<double> ctimes, stimes;

    // =====  Construction phase  =====
    /* Timers */
    double cstart, cfinish, sstart, sfinish;
    for (uint i = 0; i < cfg.ntests; i++) {
	if (B_) {
	    delete B_;
	    B_ = NULL;
	}
	cstart = cfinish = 0;
	if (cfg.prec == UH_CHEB_PREC) {
	    cstart = pclock();
	    B_ = new Prec(A, cfg);
	    cfinish = pclock();
	} else if (cfg.prec == AMG_PREC) {
	    B_ = new AMGPrec(A);
	} else if (cfg.prec == DIAG_PREC) {
	    B_ = new DiagPrec(A);
	}
	ctimes.push_back(cfinish - cstart);
	LLL_INFO("Construction time : " << ctimes.back());
    }
    PrecBase& B = *B_;

    if (cfg.prec == UH_CHEB_PREC) {
	Prec& Bcheb = static_cast<Prec&>(B);
	LOG_INFO(Bcheb);
#if 0
	Bcheb.graph_planes("grids.ps", 1, 'z', mesh);
	return 0;
#endif
    }

    double eps = 1e-6;

    Vector b(A.size(), 0.);
    /* =====  Solution phase  ===== */
    for (uint i = 0; i < cfg.ntests; i++) {
	sstart = pclock();
	if (cfg.solver == PCG_SOLVER)
	    PCGSolver(A, b, B, eps);
	else {
	    if (cfg.prec != UH_CHEB_PREC)
		THROW_EXCEPTION("Trying to call chebyshev solver for not UH preconditioner");

	    Prec& Bcheb = static_cast<Prec&>(B);
	    ChebSolver(A, Bcheb.lmin(), Bcheb.lmax(), b, Bcheb, eps);
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
