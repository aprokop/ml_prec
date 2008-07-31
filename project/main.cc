#include "config/config.h"
#include "modules/matrix/matrix.h"
#include "modules/mesh/mesh.h"

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

    Mesh mesh;
    mesh.graph_2D_plane('z', 0);

    return 0;
}
