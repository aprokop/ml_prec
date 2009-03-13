#include "main.h"
#include "config/config.h"
#include "include/time.h"
#include "modules/mesh/mesh.h"
#include "modules/prec/cheb/cheb_prec.h"
#include "modules/prec/amg/amg_prec.h"
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

    SPEMesh mesh;
    mesh.construct_matrix(A, cfg.c);

    std::cout << cfg << std::endl;

    TIME_INIT();
#if defined CHEB_PREC || defined RELX_PREC
    TIME_START();

#if defined CHEB_PREC
    Prec B(cfg.sigma, cfg.niter, A, mesh);
#else
    RelPrec B(A, cfg.niter, cfg.sigma, cfg.sigmas, cfg.mesh);
#endif

    std::cout << std::endl << TIME_INFO("Construction time") << std::endl;
    LOG_INFO(B);
    B.graph_planes("grid.ps", 1, 'z');
    exit(1);

#elif defined AMG_PREC
    AMGPrec B(A);
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
