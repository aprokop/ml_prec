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
    // Initialize logger
#ifndef NO_LOGGER
    log4cxx::PropertyConfigurator::configure("./log4cxx.properties");
#endif

    uint ncheb;
    uint nwells;
    double c, eps;
    if (set_params(argc, argv, c, eps, ncheb, nwells)) {
	LOG_INFO("Error while setting parameters, exiting...");
	std::cout << "Error while setting parameters, exiting..." << std::endl;
	return 1;
    }
    std::cout << "ncheb   = " << ncheb << std::endl;
    std::cout << "eps     = " << eps << std::endl;
    std::cout << "c       = " << c << std::endl;
    std::cout << "nwells  = " << nwells << std::endl;

    SkylineMatrix A;

    Mesh mesh(c);
    mesh.construct_matrix(A, nwells);

    TIME_INIT();
// #define CHEB_PREC
// #define AMG_PREC
#define RELX_PREC
#if   defined CHEB_PREC
    TIME_START();
    Prec B(eps, ncheb, A, mesh);
    std::cout << std::endl << TIME_INFO("Construction time") << std::endl;
    LOG_INFO(B);
    // B.graph_planes("l1.ps", 1, 'z');
    // exit(1);
#elif defined RELX_PREC
    TIME_START();
    RelPrec B(eps, A, mesh);
    std::cout << std::endl << TIME_INFO("Construction time") << std::endl;
    LOG_INFO(B);
    B.graph_planes("l1.ps", 1, 'z');
    exit(1);
#elif defined AMG_PREC
    AMGPrec B(A);
#endif

    TIME_START();
    PCG(A, Vector(A.size()), B, 1e-6);
    std::cout << std::endl << TIME_INFO("Solution time") << std::endl;

    return 0;
}
