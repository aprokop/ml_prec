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
    std::cout << "  -s|--sigmas                     Level sigmas (for cheb prec, > 1) or q (for multi-split prec, < 1)" << std::endl;
    std::cout << "  -b|--niters                     Number of iterations per level" << std::endl;
    std::cout << "  -h|--help                       Display help" << std::endl;
    std::cout << "  -c                              Value of reaction coefficient" << std::endl;
    std::cout << "  -x|--nx                         Number of points in x direction for SPE" << std::endl;
    std::cout << "  -y|--ny                         Number of points in y direction for SPE" << std::endl;
    std::cout << "  -z|--nz                         Number of points in z direction for SPE" << std::endl;
    std::cout << "  -o|--solver={cheb|pcg|simple|direct}" << std::endl;
    std::cout << "                                  Outer solver type" << std::endl;
    std::cout << "  -t|--use_tails={yes|no}         Do not use tail removing" << std::endl;
    std::cout << "  -O|--optimize-storage={yes|no}  Do not optimize storage for symmetric matrices" << std::endl;
    std::cout << "  -m|--matrix                     Matrix input file" << std::endl;
    std::cout << "  -v|--vector                     Vector input file" << std::endl;
    std::cout << "  -a|--ntests                     Number of tests to perform" << std::endl;
    std::cout << "  -u                              Construct unsymmetric matrix" << std::endl;
    std::cout << "  -S|--unsym-shift                Unsymmetric shift" << std::endl;
#ifdef HAVE_UMFPACK
    std::cout << "  -p|--prec={uh_cheb|amg|diag|gs|bgs|rbgs|sym_split|multi_split}" << std::endl;
    std::cout << "                                  Preconditioner type" << std::endl;
#else
    std::cout << "  -p|--prec={uh_cheb|amg|diag|gs|sym_split|multi_split}" << std::endl;
    std::cout << "                                  Preconditioner type" << std::endl;
#endif
    std::cout << "  -d|--dump                       Dump matrix and vector" << std::endl;
    std::cout << "     --dir                        Directory for the results (must not exist)" << std::endl;
    std::cout << "  -A|--analysis={qdropped|histogramm|q_rem_fixed_row|offdiag_ratios|1D_jacobi|col_dominance}" << std::endl;
    std::cout << "                                  Matrix analysis to perform" << std::endl;
}

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
    cfg.optimize_storage = true;
    cfg.solver           = PCG_SOLVER;
    cfg.prec             = UH_CHEB_PREC;

    cfg.dump_data       = false;
    cfg.dir	        = std::string("results/");
    cfg.analysis        = ANAL_NONE;

    static struct option long_options[] = {
	{"sigmas",		required_argument,  NULL, 's'},
	{"niters",		required_argument,  NULL, 'b'},
	{"help",		no_argument,	    NULL, 'h'},
	{"c",			required_argument,  NULL, 'c'},
	{"nx",			required_argument,  NULL, 'x'},
	{"ny",			required_argument,  NULL, 'y'},
	{"nz",			required_argument,  NULL, 'z'},
	{"solver",		required_argument,  NULL, 'o'},
	{"use_tails",		required_argument,  NULL, 't'},
	{"optimize-storage",	required_argument,  NULL, 'O'},
	{"matrix",		required_argument,  NULL, 'm'},
	{"vector",		required_argument,  NULL, 'v'},
	{"ntests",		required_argument,  NULL, 'a'},
	{"prec",		required_argument,  NULL, 'p'},
	{"dump",		no_argument,	    NULL, 'd'},
	{"unsym-shift",		required_argument,  NULL, 'S'},
	{"dir",			required_argument,  NULL, 'D'},
	{"analysis",		required_argument,  NULL, 'A'},
	{0, 0, 0, 0}
    };
    while (1) {
	int option_index = 0;
	int ch = getopt_long(argc, argv, "hb:s:O:t:c:x:y:z:to:m:v:a:p:uS:dD:A:", long_options, &option_index);

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
	    case 't': if (!strcmp(optarg, "yes"))
			  cfg.use_tails = true;
		      else if (!strcmp(optarg, "no"))
			  cfg.use_tails = false;
		      else
			  THROW_EXCEPTION("Unknown use-tails option \"" << optarg << "\"");
		      break;
	    case 'O': if (!strcmp(optarg, "yes"))
			  cfg.optimize_storage = true;
		      else if (!strcmp(optarg, "no"))
			  cfg.optimize_storage = false;
		      else
			  THROW_EXCEPTION("Unknown optimize-storage option \"" << optarg << "\"");
		      break;
	    case 'c': cfg.c = atof(optarg); break;
	    case 'x': cfg.nx = uint(atoi(optarg)); break;
	    case 'y': cfg.ny = uint(atoi(optarg)); break;
	    case 'z': cfg.nz = uint(atoi(optarg)); break;
	    case 'o': if (!strcmp(optarg, "pcg"))
			  cfg.solver = PCG_SOLVER;
		      else if (!strcmp(optarg, "cheb"))
			  cfg.solver = CHEB_SOLVER;
		      else if (!strcmp(optarg, "simple"))
			  cfg.solver = SIMPLE_SOLVER;
#ifdef HAVE_UMFPACK
		      else if (!strcmp(optarg, "direct"))
			  cfg.solver = DIRECT_SOLVER;
#endif
		      else
			  THROW_EXCEPTION("Unknown solver type \"" << optarg << "\"");
		      break;
	    case 'm': cfg.matrix = std::string(optarg); break;
	    case 'v': cfg.vector = std::string(optarg); break;
	    case 'a': cfg.ntests = uint(atoi(optarg)); break;
	    case 'p': if (!strcmp(optarg, "uh_cheb"))
			  cfg.prec = UH_CHEB_PREC;
		      else if (!strcmp(optarg, "amg"))
			  cfg.prec = AMG_PREC;
		      else if (!strcmp(optarg, "diag"))
			  cfg.prec = DIAG_PREC;
		      else if (!strcmp(optarg, "gs"))
			  cfg.prec = GS_PREC;
#ifdef HAVE_UMFPACK
		      else if (!strcmp(optarg, "bgs"))
			  cfg.prec = BGS_PREC;
		      else if (!strcmp(optarg, "rbgs"))
			  cfg.prec = RBGS_PREC;
#endif
		      else if (!strcmp(optarg, "sym_split"))
			  cfg.prec = SYM_SPLIT_PREC;
		      else if (!strcmp(optarg, "multi_split"))
			  cfg.prec = MULTI_SPLIT_PREC;
		      else
			  THROW_EXCEPTION("Unknown solver type \"" << optarg << "\"");
		      break;
	    case 'u': cfg.unsym_matrix = true; break;
	    case 'S': cfg.unsym_shift = atof(optarg); break;
	    case 'd': cfg.dump_data = true; break;
	    case 'D': cfg.dir = std::string(optarg); break;
	    case 'A': if (!strcmp(optarg, "qdropped"))
			  cfg.analysis = ANAL_QDROPPED;
		      else if (!strcmp(optarg, "none"))
			  cfg.analysis = ANAL_NONE;
		      else if (!strcmp(optarg, "histogramm"))
			  cfg.analysis = ANAL_HISTOGRAMM;
		      else if (!strcmp(optarg, "q_rem_fixed_row"))
			  cfg.analysis = ANAL_Q_REM_FIXED_ROW;
		      else if (!strcmp(optarg, "offdiag_ratios"))
			  cfg.analysis = ANAL_OFFDIAGONAL_RATIOS;
		      else if (!strcmp(optarg, "1D_jacobi"))
			  cfg.analysis = ANAL_1D_JACOBI;
		      else if (!strcmp(optarg, "col_dominance"))
			  cfg.analysis = ANAL_COL_DOMINANCE;
		      else
			  THROW_EXCEPTION("Unknown analysis type \"" << optarg << "\"");
		      break;
	    case '?':
	    default:
		      abort();
	}
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
    if (cfg.prec == MULTI_SPLIT_PREC)
	for (uint i = 0; i < cfg.sigmas.size(); i++)
	    if (cfg.sigmas[i] >= 1)
		THROW_EXCEPTION("All sigmas must be < 1");

    if (cfg.unsym_matrix == true)
	cfg.solver = SIMPLE_SOLVER;

    /* Create directory for results */
    // if (mkdir("res", 0700) == -1) {
	// THROW_EXCEPTION("Cannot create directory \"" << cfg.dir << "\" (" << strerror(errno) << ")");
	// return 1;
    // }

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
	if (!cfg.unsym_matrix)
	    mesh.construct_matrix(A, cfg.c);
	else
	    mesh.construct_matrix_unsym(A, cfg.c, cfg.unsym_shift);
    } else {
	/*
	 * Read the matrix.
	 * If matrix is written in CSR format we'll need to convert it to Skyline (transform = true)
	 */
	bool transform = true;

	DumpType type = ((strstr(cfg.matrix.c_str(), ".crs") ? true : false)) ? ASCII : BINARY;
	A.load(cfg.matrix, transform, type);
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
#else
    /* Dump matrix in the BINARY format */
    dump("matrix_hypre.dat", A, BINARY);
#endif
    dump("vector_hypre.dat.00000", b, HYPRE);
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

    if (cfg.analysis != ANAL_NONE) {
	os << "Analysis         : ";
	switch(cfg.analysis) {
	    case ANAL_HISTOGRAMM         : os << "histogramm"; break;
	    case ANAL_QDROPPED           : os << "qdropped"; break;
	    case ANAL_Q_REM_FIXED_ROW    : os << "q_rem_fixed_row"; break;
	    case ANAL_OFFDIAGONAL_RATIOS : os << "offdiag_ratios"; break;
	    case ANAL_1D_JACOBI          : os << "1D_jacobi";
	    case ANAL_COL_DOMINANCE	 : os << "col_dominance"; break;
	    case ANAL_NONE		 : break;
	}
	os << std::endl;

	return os;
    }

    /* Run parameters */
    os << "Number of tests  : " << cfg.ntests << std::endl;
    os << "Use tails        : " << (cfg.use_tails ? "true" : "false") << std::endl;
    os << "Optimize storage : " << (cfg.optimize_storage ? "true" : "false") << std::endl;

    os << "Solver           : ";
    switch(cfg.solver) {
	case PCG_SOLVER    : os << "pcg"; break;
	case CHEB_SOLVER   : os << "cheb"; break;
	case SIMPLE_SOLVER : os << "simple"; break;
#ifdef HAVE_UMFPACK
	case DIRECT_SOLVER : os << "direct"; break;
#endif
	default		: THROW_EXCEPTION("Unknown SolverType");
    }
    os << std::endl;

    os << "Preconditioner   : ";
    switch(cfg.prec) {
	case UH_CHEB_PREC     : os << "uh_cheb"; break;
	case AMG_PREC         : os << "amg"; break;
	case DIAG_PREC        : os << "diag"; break;
	case GS_PREC          : os << "gs"; break;
#ifdef HAVE_UMFPACK
	case BGS_PREC         : os << "bgs"; break;
	case RBGS_PREC        : os << "rbgs"; break;
#endif
	case SYM_SPLIT_PREC   : os << "sym_split"; break;
	case MULTI_SPLIT_PREC : os << "multi_split"; break;
	default           : THROW_EXCEPTION("Unknown PrecType");
    }
    os << std::endl;

    return os;
}

void GlobalStats::dump(const std::string& dir) const {
    std::ofstream ofs_const((dir + "t_const").c_str()); ofs_const << std::fixed << std::setprecision(3) << t_const; ofs_const.close();
    std::ofstream ofs_prec((dir + "t_prec").c_str());   ofs_prec  << std::fixed << std::setprecision(3) << t_prec;  ofs_prec.close();
    std::ofstream ofs_sol((dir + "t_sol").c_str());     ofs_sol   << std::fixed << std::setprecision(3) << t_sol;   ofs_sol.close();
    std::ofstream ofs_total((dir + "t_total").c_str()); ofs_total << std::fixed << std::setprecision(3) << t_total; ofs_total.close();
    std::ofstream ofs_niter((dir + "niter").c_str());   ofs_niter << std::fixed << std::setprecision(3) << niter;   ofs_niter.close();
}
