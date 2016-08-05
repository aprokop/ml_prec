#include "config.h"
#include "include/logger.h"
#include "include/time.h"
#include "include/tools.h"
#include "main.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iterator>
#include <numeric>
#include <fstream>

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

DEFINE_LOGGER("Main");

static void usage() {
    std::cout << "Usage: ./spe_prec [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -a|--ntests                     Number of tests to perform" << std::endl;
    std::cout << "  -A|--analysis={qdropped|histogramm|q_rem_fixed_row|offdiag_ratios|1D_jacobi|col_dominance|unsym_convergence}" << std::endl;
    std::cout << "                                  Matrix analysis to perform" << std::endl;
    std::cout << "  -b|--niters                     Number of iterations per level" << std::endl;
    std::cout << "  -c                              Value of reaction coefficient" << std::endl;
    std::cout << "  -C|--prec-config-file           Config file for the preconditioner" << std::endl;
    std::cout << "  -d|--dump                       Dump matrix and vector" << std::endl;
    std::cout << "     --dir                        Directory for the results (must not exist)" << std::endl;
    std::cout << "  -h|--help                       Display help" << std::endl;
    std::cout << "  -m|--matrix                     Matrix input file" << std::endl;
#ifdef HAVE_UMFPACK
    std::cout << "  -M|--max-levels                 Maximum number of levels for constructed preconditioner" << std::endl;
    std::cout << "  -N|--coarse-n                   Use director solver for coarse systems of order <= N" << std::endl;
#endif
    std::cout << "  -o|--solver={cheb|pcg|simple|gmres|direct}" << std::endl;
    std::cout << "                                  Outer solver type" << std::endl;
    std::cout << "  -O|--optimize-storage={yes|no}  Do not optimize storage for symmetric matrices" << std::endl;
#ifdef HAVE_UMFPACK
    std::cout << "  -p|--prec={uh_cheb|amg|comp|diag|gs|id|bgs|rbgs|sym_split|multi_split}" << std::endl;
    std::cout << "                                  Preconditioner type" << std::endl;
#else
    std::cout << "  -p|--prec={uh_cheb|amg|comp|diag|gs|id|sym_split|multi_split}" << std::endl;
    std::cout << "                                  Preconditioner type" << std::endl;
#endif
    std::cout << "  -s|--sigmas                     Level sigmas (for cheb prec, > 1) or q (for multi-split prec, < 1)" << std::endl;
    std::cout << "  -S|--unsym-shift                Unsymmetric shift" << std::endl;
    std::cout << "  -t|--use_tails={yes|no}         Do not use tail removing" << std::endl;
    std::cout << "  -T|--transform={none|IL|IU|ILU}" << std::endl;
    std::cout << "                                  Transformation to perform" << std::endl;
    std::cout << "  -u                              Construct unsymmetric matrix" << std::endl;
    std::cout << "  -v|--vector                     Vector input file" << std::endl;
    std::cout << "  -x|--nx                         Number of points in x direction for SPE" << std::endl;
    std::cout << "  -y|--ny                         Number of points in y direction for SPE" << std::endl;
    std::cout << "  -z|--nz                         Number of points in z direction for SPE" << std::endl;
}

static const uint DEFAULT_MAX_LEVELS = 30;
int set_params(int argc, char * argv[], Config& cfg) {
    cfg.niters.push_back(2);
    cfg.sigmas.push_back(3);

    cfg.c                = 1.;
    cfg.nx               = 60;
    cfg.ny               = 220;
    cfg.nz               = 85;

    cfg.unsym_matrix	 = false;
    cfg.unsym_shift	 = 0.1;

    cfg.ntests           = 1;
    cfg.use_tails        = true;
    cfg.coarse_n	 = 0;
    cfg.max_levels	 = DEFAULT_MAX_LEVELS;
    cfg.optimize_storage = true;
    cfg.solver           = PCG_SOLVER;
    cfg.prec             = UH_CHEB_PREC;

    cfg.dump_data        = false;
    cfg.dir	         = std::string("results/");
    cfg.analysis         = ANAL_NONE;
    cfg.transform	 = TRANS_NONE;

    static struct option long_options[] = {
        {"ntests",		required_argument,  NULL, 'a'},
        {"analysis",		required_argument,  NULL, 'A'},
        {"niters",		required_argument,  NULL, 'b'},
        {"c",			required_argument,  NULL, 'c'},
        {"prec-config-file",	required_argument,  NULL, 'C'},
        {"dump",		no_argument,	    NULL, 'd'},
        {"dir",			required_argument,  NULL, 'D'},
        {"help",		no_argument,	    NULL, 'h'},
        {"matrix",		required_argument,  NULL, 'm'},
        {"max-levels",		required_argument,  NULL, 'M'},
        {"coarse-n",		required_argument,  NULL, 'N'},
        {"optimize-storage",	required_argument,  NULL, 'O'},
        {"solver",		required_argument,  NULL, 'o'},
        {"prec",		required_argument,  NULL, 'p'},
        {"sigmas",		required_argument,  NULL, 's'},
        {"unsym-shift",		required_argument,  NULL, 'S'},
        {"use_tails",		required_argument,  NULL, 't'},
        {"transform",		required_argument,  NULL, 'T'},
        {"vector",		required_argument,  NULL, 'v'},
        {"nx",			required_argument,  NULL, 'x'},
        {"ny",			required_argument,  NULL, 'y'},
        {"nz",			required_argument,  NULL, 'z'},
        {0, 0, 0, 0}
    };

#define CHECK_AND_SET(str, param, value) if (!strcmp(optarg, str)) param = value
    while (1) {
        int option_index = 0;
        int ch = getopt_long(argc, argv, "s:O:a:A:b:c:C:dD:hm:M:N:o:p:S:t:T:uv:x:y:z:", long_options, &option_index);

        if (ch == -1)
            break;

        if (ch == 0)
            continue;

        switch (ch) {
            case 'h': usage();
                      exit(0);
            case 's': {
                          cfg.sigmas.clear();
                          std::istringstream os(optarg);
                          std::copy(std::istream_iterator<double>(os), std::istream_iterator<double>(),
                                    std::back_inserter(cfg.sigmas));
                      }
                      break;
            case 'b': {
                          cfg.niters.clear();
                          std::istringstream os(optarg);
                          std::copy(std::istream_iterator<double>(os), std::istream_iterator<double>(),
                                    std::back_inserter(cfg.niters));
                      }
                      break;
            case 't': CHECK_AND_SET("yes", cfg.use_tails, true);
                      else CHECK_AND_SET("no", cfg.use_tails, false);
                      else THROW_EXCEPTION("Unknown use-tails option \"" << optarg << "\"");
                      break;
            case 'O': CHECK_AND_SET("yes", cfg.optimize_storage, true);
                      else CHECK_AND_SET("no", cfg.optimize_storage, false);
                      else THROW_EXCEPTION("Unknown optimize-storage option \"" << optarg << "\"");
                      break;
            case 'c': cfg.c = atof(optarg); break;
            case 'x': cfg.nx = uint(atoi(optarg)); break;
            case 'y': cfg.ny = uint(atoi(optarg)); break;
            case 'z': cfg.nz = uint(atoi(optarg)); break;
            case 'M': cfg.max_levels = uint(atoi(optarg)); break;
            case 'N': cfg.coarse_n = uint(atoi(optarg)); break;
            case 'o': CHECK_AND_SET("pcg", cfg.solver, PCG_SOLVER);
                      else CHECK_AND_SET("cheb", cfg.solver, CHEB_SOLVER);
                      else CHECK_AND_SET("simple", cfg.solver, SIMPLE_SOLVER);
                      else CHECK_AND_SET("gmres", cfg.solver, GMRES_SOLVER);
#ifdef HAVE_UMFPACK
                      else CHECK_AND_SET("direct", cfg.solver, DIRECT_SOLVER);
#endif
                      else THROW_EXCEPTION("Unknown solver type \"" << optarg << "\"");
                      break;
            case 'm': cfg.matrix = std::string(optarg); break;
            case 'v': cfg.vector = std::string(optarg); break;
            case 'a': cfg.ntests = uint(atoi(optarg)); break;
            case 'p': CHECK_AND_SET("uh_cheb", cfg.prec, UH_CHEB_PREC);
                      else CHECK_AND_SET("amg", cfg.prec, AMG_PREC);
                      else CHECK_AND_SET("comp", cfg.prec, COMP_PREC);
                      else CHECK_AND_SET("diag", cfg.prec, DIAG_PREC);
                      else CHECK_AND_SET("gs", cfg.prec, GS_PREC);
                      else CHECK_AND_SET("id", cfg.prec, ID_PREC);
#ifdef HAVE_UMFPACK
                      else CHECK_AND_SET("bgs", cfg.prec, BGS_PREC);
                      else CHECK_AND_SET("rbgs", cfg.prec, RBGS_PREC);
#endif
                      else CHECK_AND_SET("sym_split", cfg.prec, SYM_SPLIT_PREC);
                      else CHECK_AND_SET("multi_split", cfg.prec, MULTI_SPLIT_PREC);
                      else THROW_EXCEPTION("Unknown prec type \"" << optarg << "\"");
                      break;
            case 'u': cfg.unsym_matrix = true; break;
            case 'S': cfg.unsym_shift = atof(optarg); break;
            case 'd': cfg.dump_data = true; break;
            case 'C': cfg.prec_conf_file = std::string(optarg); break;
            case 'D': cfg.dir = std::string(optarg); break;
            case 'A': CHECK_AND_SET("none", cfg.analysis, ANAL_NONE);
                      else CHECK_AND_SET("qdropped", cfg.analysis, ANAL_QDROPPED);
                      else CHECK_AND_SET("histogramm", cfg.analysis, ANAL_HISTOGRAMM);
                      else CHECK_AND_SET("q_rem_fixed_row", cfg.analysis, ANAL_Q_REM_FIXED_ROW);
                      else CHECK_AND_SET("offdiag_ratios", cfg.analysis, ANAL_OFFDIAGONAL_RATIOS);
                      else CHECK_AND_SET("1D_Jacobi", cfg.analysis, ANAL_1D_JACOBI);
                      else CHECK_AND_SET("col_dominance", cfg.analysis, ANAL_COL_DOMINANCE);
                      else CHECK_AND_SET("unsym_convergence", cfg.analysis, ANAL_2LEVEL_CONVERGENCE);
                      else THROW_EXCEPTION("Unknown analysis type \"" << optarg << "\"");
                      break;
            case 'T': CHECK_AND_SET("none", cfg.transform, TRANS_NONE);
                      else CHECK_AND_SET("IL", cfg.transform, TRANS_IL);
                      else CHECK_AND_SET("IU", cfg.transform, TRANS_IU);
                      else CHECK_AND_SET("ILU", cfg.transform, TRANS_ILU);
                      else THROW_EXCEPTION("Unknown transformation type \"" << optarg << "\"");
                      break;
            case '?':
            default:
                      abort();
        }
    }
#undef CHECK_AND_SET

#ifndef HAVE_UMFPACK
    if (cfg.coarse_n) {
        cfg.coarse_n = 0;
        LLL_WARN("No direct solver available, ignoring -N argument");
    }
    if (cfg.max_levels) {
        cfg.max_levels = 0;
        LLL_WARN("No direct solver available, ignoring -M argument");
    }
#endif

    if (cfg.analysis == ANAL_2LEVEL_CONVERGENCE && cfg.prec != MULTI_SPLIT_PREC) {
        LLL_WARN("Changing preconditioner to \"multi_split\" for this type of analysis");
        cfg.prec = MULTI_SPLIT_PREC;
    }

    /* Check inconsistencies in parameters */
    if (cfg.prec == UH_CHEB_PREC) {
        if (cfg.niters.size() == 1 && cfg.sigmas.size() == 1 &&
            ((cfg.niters[0] == 2 && cfg.sigmas[0] >= 4) ||
             (cfg.niters[0] == 3 && cfg.sigmas[0] >= 9))) {
            LOG_WARN("Possibly wrong relation between niter and sigma");
        }
    }
    if (cfg.solver == CHEB_SOLVER && cfg.prec != UH_CHEB_PREC)
        THROW_EXCEPTION("Trying to call Chebyshev solver for not UH Cheb preconditioner");

    if (cfg.prec == UH_CHEB_PREC || cfg.prec == SYM_SPLIT_PREC)
        for (uint i = 0; i < cfg.sigmas.size(); i++)
            if (cfg.sigmas[i] <= 1)
                THROW_EXCEPTION("All sigmas must be > 1");
    if (cfg.prec == MULTI_SPLIT_PREC && cfg.prec_conf_file.empty())
        for (uint i = 0; i < cfg.sigmas.size(); i++)
            if (cfg.sigmas[i] >= 1)
                THROW_EXCEPTION("All sigmas must be < 1");

    if (cfg.unsym_matrix == true &&
        (cfg.solver != SIMPLE_SOLVER && cfg.solver != GMRES_SOLVER))
        cfg.solver = SIMPLE_SOLVER;

    return 0;
}

double avg_time(const std::vector<double>& times) {
    const uint n = times.size();
    std::vector<double> ts = times;
    ASSERT(n, "zero length vector");

    std::sort(ts.begin(), ts.end());
    if (n <= 3)
        return ts[0];

    /* We ignore two largest ts and 1 smallest and average the rest */
    return std::accumulate(ts.begin() + 1, ts.end() - 2, 0.) / (n-3);
}

void construct_matrix(const Config& cfg, const SPEMesh& mesh, SkylineMatrix& A) {
    if (cfg.matrix.empty()) {
        /* Construct the matrix */
        if (!cfg.unsym_matrix)	mesh.construct_matrix(A, cfg.c);
        else			mesh.construct_matrix_unsym(A, cfg.c, cfg.unsym_shift);
    } else {
        /*
         * Read the matrix.
         * If matrix is written in CSR format we'll need to convert it to Skyline (transform = true)
         */
        bool transform = true;

        DumpType type = BINARY;
        if (strstr(cfg.matrix.c_str(), ".crs")) type = ASCII;
        if (strstr(cfg.matrix.c_str(), ".mm"))  type = MATRIX_MARKET;
        A.load(cfg.matrix, type, transform);
    }
}

void construct_vector(const Config& cfg, Vector& b) {
    if (!cfg.vector.empty()) {
        DumpType type = ((strstr(cfg.vector.c_str(), ".crs") ? true : false)) ? ASCII : BINARY;
        load(b, cfg.vector, type);
    }
}

void dump_data(const SkylineMatrix& A, const Vector& b) {
#if 0
    /* Dump matrix in the HYPRE format for further running of HYPRE */
    dump("matrix_hypre.dat.00000", A, HYPRE);
    dump("vector_hypre.dat.00000", b, HYPRE);
#else
    /* Dump matrix in the BINARY format */
    dump("matrix.dat", A, BINARY);
    dump("vector.dat", b, BINARY);
#endif
}

std::ostream& operator<<(std::ostream& os, const Config& cfg) {
    os << "Level iterations : ";
    for (uint i = 0; i < cfg.niters.size(); i++)
        os << cfg.niters[i] << " ";
    os << std::endl;
    os << "Level sigmas     : ";
    for (uint i = 0; i < cfg.sigmas.size(); i++)
        os << cfg.sigmas[i] << " ";
    os << std::endl;

    /* Mesh parameters */
    if (cfg.matrix.empty()) {
        os << "c                : " << cfg.c << std::endl;
        os << "Geometry         : " << cfg.nx << " x " << cfg.ny << " x " << cfg.nz << std::endl;
        os << "Unsymmetric      : " << (cfg.unsym_matrix ? "true" : "false") << std::endl;
        if (cfg.unsym_matrix)
            os << "Unsymmetrix shift: " << cfg.unsym_shift << std::endl;
    } else {
        os << "Matrix file      : " << cfg.matrix << std::endl;
    }
    os << "Dump data        : " << (cfg.dump_data ? "true" : "false") << std::endl;

#define CASE_PRINT(value, str) case value: os << str; break
    if (cfg.analysis != ANAL_NONE) {
        os << "Analysis         : ";
        switch(cfg.analysis) {
            CASE_PRINT(ANAL_HISTOGRAMM, "histogramm");
            CASE_PRINT(ANAL_QDROPPED, "qdropped");
            CASE_PRINT(ANAL_Q_REM_FIXED_ROW, "q_rem_fixed_row");
            CASE_PRINT(ANAL_OFFDIAGONAL_RATIOS, "offdiag_ratios");
            CASE_PRINT(ANAL_1D_JACOBI, "1D_jacobi");
            CASE_PRINT(ANAL_COL_DOMINANCE, "col_dominance");
            CASE_PRINT(ANAL_2LEVEL_CONVERGENCE, "unsym_convergence");
            case ANAL_NONE: break;
        }
        os << std::endl;

        return os;
    }

    if (cfg.transform != TRANS_NONE) {
        os << "Transformation   : ";
        switch(cfg.transform) {
            CASE_PRINT(TRANS_IL, "IL");
            CASE_PRINT(TRANS_IU, "IU");
            CASE_PRINT(TRANS_ILU, "ILU");
            case TRANS_NONE: break;
        }
        os << std::endl;

        return os;
    }

    /* Run parameters */
    os << "Number of tests  : " << cfg.ntests << std::endl;
    os << "Use tails        : " << (cfg.use_tails ? "true" : "false") << std::endl;
    os << "Coarse n         : " << cfg.coarse_n << std::endl;
    if (cfg.max_levels != DEFAULT_MAX_LEVELS &&
        (cfg.prec == UH_CHEB_PREC || cfg.prec == MULTI_SPLIT_PREC))
        os << "Max levels       : " << cfg.max_levels << std::endl;
    os << "Optimize storage : " << (cfg.optimize_storage ? "true" : "false") << std::endl;

    os << "Solver           : ";
    switch(cfg.solver) {
        CASE_PRINT(PCG_SOLVER, "pcg");
        CASE_PRINT(CHEB_SOLVER, "cheb");
        CASE_PRINT(SIMPLE_SOLVER, "simple");
        CASE_PRINT(GMRES_SOLVER, "gmres");
#ifdef HAVE_UMFPACK
        CASE_PRINT(DIRECT_SOLVER, "direct");
#endif
        default		: THROW_EXCEPTION("Unknown SolverType");
    }
    os << std::endl;

    os << "Preconditioner   : ";
    switch(cfg.prec) {
        CASE_PRINT(AMG_PREC, "amg");
        CASE_PRINT(COMP_PREC, "composite");
        CASE_PRINT(DIAG_PREC, "diag");
        CASE_PRINT(GS_PREC, "gs");
        CASE_PRINT(ID_PREC, "id");
#ifdef HAVE_UMFPACK
        CASE_PRINT(BGS_PREC, "bgs");
        CASE_PRINT(RBGS_PREC, "rbgs");
#endif
        CASE_PRINT(MULTI_SPLIT_PREC, "multi_split");
        CASE_PRINT(SYM_SPLIT_PREC, "sym_split");
        CASE_PRINT(UH_CHEB_PREC, "uh_cheb");
        default: THROW_EXCEPTION("Unknown PrecType");
    }
    os << std::endl;

#undef CASE_PRINT

    return os;
}

void GlobalStats::dump(const std::string& dir) const {
    std::ofstream ofs_const((dir + "t_const").c_str()); ofs_const << std::fixed << std::setprecision(3) << t_const; ofs_const.close();
    std::ofstream ofs_prec((dir + "t_prec").c_str());   ofs_prec  << std::fixed << std::setprecision(3) << t_prec;  ofs_prec.close();
    std::ofstream ofs_sol((dir + "t_sol").c_str());     ofs_sol   << std::fixed << std::setprecision(3) << t_sol;   ofs_sol.close();
    std::ofstream ofs_total((dir + "t_total").c_str()); ofs_total << std::fixed << std::setprecision(3) << t_total; ofs_total.close();
    std::ofstream ofs_niter((dir + "niter").c_str());   ofs_niter << std::fixed << std::setprecision(3) << niter;   ofs_niter.close();
}
