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

// #define CHEB_PREC
// #define AMG_PREC
#define RELX_PREC

int main (int argc, char * argv[]) {
#ifndef NO_LOGGER
    // Initialize logger
    log4cxx::PropertyConfigurator::configure("./log4cxx.properties");
#endif

    uint niter;
    uint nwells;
    double c, sigma, gamma;
    std::vector<double> sigmas;
    if (set_params(argc, argv, c, sigma, sigmas, niter, nwells)) {
	LOG_INFO("Error while setting parameters, exiting...");
	std::cout << "Error while setting parameters, exiting..." << std::endl;
	return 1;
    }
    std::cout << "c       = " << c << std::endl;
#if   defined CHEB_PREC
    std::cout << "ncheb   = " << niter << std::endl;
    std::cout << "sigma   = " << sigma << std::endl;
#elif defined RELX_PREC
    std::cout << "niter   = " << niter << std::endl;
    std::cout << "gamma   = " << sigma << std::endl;
    std::cout << "sigmas  = ";
    for (uint i = 0; i < sigmas.size(); i++)
	std::cout << sigmas[i] << " ";
    std::cout << std::endl;
#endif
    // std::cout << "nwells  = " << nwells << std::endl;

    SkylineMatrix A;

    SPEMesh mesh;
    mesh.construct_matrix(A, c);

    TIME_INIT();
#if defined CHEB_PREC || defined RELX_PREC
    TIME_START();

#if defined CHEB_PREC
    Prec B(sigma, niter, A, mesh);
#else
    RelPrec B(A, niter, sigma, sigmas, mesh);
#endif

    std::cout << std::endl << TIME_INFO("Construction time") << std::endl;
    LOG_INFO(B);
    // B.graph_planes("grid.ps", 1, 'z');
    // exit(1);

#elif defined AMG_PREC
    AMGPrec B(A);
#endif

    TIME_START();
    PCG(A, Vector(A.size()), B, 1e-6);
    std::cout << std::endl << TIME_INFO("Solution time") << std::endl;

    return 0;
}
