#include "main.h"
#include "config/config.h"
#include "include/time.h"
#include "modules/matrix/matrix.h"
#include "modules/mesh/mesh.h"
#include "modules/prec/amg/amg_prec.h"
#include "modules/prec/diag/diag.h"
#include "modules/prec/cheb/cheb_prec.h"
#include "modules/prec/relax/rel_prec.h"
#include "modules/prec/sym/sym.h"
#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/solvers/solvers.h"

/* Logger header files */
#include "include/logger.h"
#ifndef NO_LOGGER
#include <log4cxx/propertyconfigurator.h>
#endif

DEFINE_LOGGER("Main");

int main (int argc, char * argv[]) {
#ifndef NO_LOGGER
    /* Initialize logger */
    log4cxx::PropertyConfigurator::configure("./log4cxx.properties");
#endif

    /* Set problem parameters */
    Config cfg;
    if (set_params(argc, argv, cfg)) {
        LLL_INFO("Error while setting parameters, exiting...");
        return 1;
    }

    LLL_DEBUG("Config parameters: " << cfg);


    SPEMesh mesh(cfg.nx, cfg.ny, cfg.nz);

    SkylineMatrix A;
    Vector b(A.size(), 0.);
    construct_matrix(cfg, mesh, A);
    construct_vector(cfg, b);

    const double IV = 1000; /* Initial value */

#define INDEX(i,j,k) ((k)*220*60 + (j)*60 + (i))
    /* Set some initial values to IV to examine the spread */
    b[INDEX(13,65,55)]  =
            b[INDEX(48,77,55)]  =
            b[INDEX(48,77,10)]  =
            b[INDEX(17,150,55)] = cfg.c*IV;
#undef INDEX

    DiagPrec B(A);

    /* We need to solve the system exactly, so more precise eps */
    double eps = 1e-13;

    Vector x(A.size());
    if (!cfg.unsym_matrix) {
        if (cfg.solver == PCG_SOLVER)
            PCGSolver(A, b, B, x, eps);
        else {
            /* ChebSolver requires lower and upper bounds on eigenvalues. Thus
             * it is applicable only to Chebyshev preconditioner */
            Prec& Bcheb = static_cast<Prec&>(B);
            ChebSolver(A, Bcheb.lmin(), Bcheb.lmax(), b, Bcheb, x, eps);
        }

        /* Dump solution data for fading for further examination (via Python scripts) */
        std::ofstream os("vector.bin", std::ofstream::binary);
        uint n = x.size();
        for (uint i = 0; i < n; i++) {
            double d = x[i];
            if (d < 0) {
                LOG_WARN("x[" << i << "] = " << x[i]);
                x[i] = 0;
            }
            // d /= IV;
            d = log(x[i]/IV + 1e-15)/log(10);
            os.write(reinterpret_cast<const char*>(&d), sizeof(double));
        }
        os.close();

        return 0;
    }
