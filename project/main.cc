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

    // logger.logf(log4cxx::Level::DEBUG, "Some forma: %d %%\n", 15);

    Config cfg;
    if (set_params(argc, argv, cfg)) {
	LOG_INFO("Error while setting parameters, exiting...");
	std::cout << "Error while setting parameters, exiting..." << std::endl;
	return 1;
    }

    SkylineMatrix A;

    SPEMesh mesh;
    mesh.construct_matrix(A, cfg.c);
    A.stat(true);

    std::cout << cfg << std::endl;
    LOG_DEBUG("Config parameters: " << cfg);

    TIME_INIT();

    PrecBase * B_ = 0;
    if (cfg.prec == UH_CHEB_PREC) {
	TIME_START();
	B_ = new Prec(A, cfg);
	std::cout << std::endl << TIME_INFO("Construction time") << std::endl;

	Prec& Bcheb = static_cast<Prec&>(*B_);
	LOG_INFO(Bcheb);
#if 1
	Bcheb.graph_planes("grids.ps", 0, 'z', mesh);
	return 0;
#endif
    } else if (cfg.prec == AMG_PREC) {
	B_ = new AMGPrec(A);
    } else if (cfg.prec == DIAG_PREC) {
	B_ = new DiagPrec(A);
    }
    PrecBase& B = *B_;

    TIME_START();
    if (cfg.solver == PCG_SOLVER)
	PCGSolver(A, Vector(A.size()), B, 1e-6);
    else {
	if (cfg.prec != UH_CHEB_PREC)
	    THROW_EXCEPTION("Trying to call chebyshev solver for not UH preconditioner");

	Prec& Bcheb = static_cast<Prec&>(B);
	ChebSolver(A, Bcheb.lmin(), Bcheb.lmax(), Vector(A.size()), Bcheb, 1e-6);
    }
    std::cout << std::endl << TIME_INFO("Solution time") << std::endl;

    delete B_;

    return 0;
}
