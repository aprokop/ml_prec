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

// logger
#include "include/logger.h"
#ifndef NO_LOGGER
#include <log4cxx/propertyconfigurator.h>
#endif

DEFINE_LOGGER("Main");

#define TEST_FADING

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

    std::cout << cfg << std::endl;
    LOG_DEBUG("Config parameters: " << cfg);

    SkylineMatrix A;

    SPEMesh mesh(cfg.nx, cfg.ny, cfg.nz);
    if (cfg.matrix.empty()) {
	/* Construct the matrix */
	if (!cfg.unsym_matrix)
	    mesh.construct_matrix(A, cfg.c);
	else
	    mesh.construct_matrix_unsym(A, cfg.c, cfg.unsym_shift);
    } else {
	/*
	 * Read the matrix.
	 * If matrix is writtent in CSR format we'll need to convert it to Skyline (transform = true)
	 */
	bool transform = true;
	A.load(cfg.matrix, transform);
    }

    Vector b(A.size(), 0.);
    if (!cfg.vector.empty()) {
	/* By default, we load vector in ASCII mode, cause binary is not implemented yet */
	load(b, cfg.vector, ASCII);
    }

#ifdef TEST_FADING
    const double IV = 1000; /* Initial value */

#define INDEX(i,j,k) ((k)*220*60 + (j)*60 + (i))
    /* Place several probes to some points */
    b[INDEX(13,65,55)]  =
    b[INDEX(48,77,55)]  =
    b[INDEX(48,77,10)]  =
    b[INDEX(17,150,55)] = cfg.c*IV;
#undef INDEX
#endif

    if (cfg.dump_data) {
#if 0
	dump("matrix_hypre.dat.00000", A, HYPRE);
#else
	dump("matrix_hypre.dat", A, BINARY);
#endif
	dump("vector_hypre.dat.00000", b, HYPRE);
	exit(0);
    }

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
	cstart = pclock();
	switch (cfg.prec) {
	    case UH_CHEB_PREC     : B_ = new Prec(A, cfg);		break;
	    case AMG_PREC         : B_ = new AMGPrec(A);		break;
	    case DIAG_PREC        : B_ = new DiagPrec(A);		break;
	    case SYM_SPLIT_PREC   : B_ = new SymPrec(A, cfg);		break;
	    case MULTI_SPLIT_PREC : B_ = new MultiSplitPrec(A, cfg);	break;
	}
	cfinish = pclock();
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
    } else if (cfg.prec == MULTI_SPLIT_PREC) {
	MultiSplitPrec& Bms = static_cast<MultiSplitPrec&>(B);
	LOG_INFO(Bms);
#if 0
	Bms.graph_planes("grids.ps", 1, 'z', mesh);
	return 0;
#endif
    }

    double eps = 1e-6;
#ifdef TEST_FADING
    eps = 1e-13;
#endif

    Vector x(A.size());
    /* =====  Solution phase  ===== */
    for (uint i = 0; i < cfg.ntests; i++) {
	sstart = pclock();
	if (!cfg.unsym_matrix) {
	    if (cfg.solver == PCG_SOLVER)
		PCGSolver(A, b, B, x, eps);
	    else {
		Prec& Bcheb = static_cast<Prec&>(B);
		ChebSolver(A, Bcheb.lmin(), Bcheb.lmax(), b, Bcheb, x, eps);
	    }
	} else {
	    SimpleSolver(A, b, B, x, eps);
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

#ifdef TEST_FADING
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
#endif


    return 0;
}
