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
#if defined CHEB_PREC || defined RELX_PREC
    TIME_START();

#if defined CHEB_PREC
    Prec B(A, cfg);
#else
    RelPrec B(A, cfg.niter, cfg.sigma, cfg.sigmas, cfg.mesh);
#endif

    std::cout << std::endl << TIME_INFO("Construction time") << std::endl;
    LOG_INFO(B);

#elif defined AMG_PREC
    AMGPrec B(A);
#elif defined DIAG_PREC
    DiagPrec B(A);
#endif

    TIME_START();
#if 1
    PCGSolver(A, Vector(A.size()), B, 1e-6);
#else
#ifndef CHEB_PREC
#  error "CHEB_PREC is not defined"
#endif
    ChebSolver(A, B.lmin(), B.lmax(), Vector(A.size()), B, 1e-6);
#endif
    std::cout << std::endl << TIME_INFO("Solution time") << std::endl;

    return 0;
}
