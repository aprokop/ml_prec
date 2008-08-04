#include "main.h"
#include "config/config.h"
#include "include/time.h"
#include "modules/mesh/mesh.h"
#include "modules/prec/prec.h"
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

    double c = 1.0;

    uint nlevels, ncheb;
    double eps;
    if (set_params(argc, argv, nlevels, eps, ncheb)) {
	LOG_INFO("Error while setting parameters, exiting...");
	std::cout << "Error while setting parameters, exiting..." << std::endl;
	return 1;
    }
    std::cout << "nlevels = " << nlevels << std::endl;
    std::cout << "eps = " << eps << std::endl;
    std::cout << "ncheb = " << ncheb << std::endl;
    std::cout << "c = " << c << std::endl;

    Mesh mesh(c);
    const CSRMatrix& A = mesh.get_matrix();
    // mesh.graph_xy_planes();
    // mesh.graph_z_lines();

    TIME_INIT();
#if 0
    TIME_START();
    Prec B(nlevels, eps, ncheb, c, A);
    std::cout << std::endl << TIME_INFO("Construction time: ") << std::endl;
    LOG_INFO(B);
#else
    AMGPrec B(A);
#endif

    PCG(A, Vector(A.size()), B, 1e-6);

    return 0;
}
