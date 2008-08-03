#include "config/config.h"
#include "main.h"
#include "modules/mesh/mesh.h"
#include "modules/prec/prec.h"

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
    // mesh.graph_xy_planes();
    // mesh.graph_z_lines();

    Prec B(nlevels, eps, ncheb, c, mesh.get_matrix());
    LOG_INFO(B);

    return 0;
}
