#include "config/config.h"
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

DEFINE_LOGGER("Main");

static void usage() {
    std::cout << "Usage: ./spe_prec [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -s|--sigmas                     Level sigmas" << std::endl;
    std::cout << "  -b|--niters                     Number of iterations per level" << std::endl;
    std::cout << "  -h|--help                       Display help" << std::endl;
    std::cout << "  -c                              Value of reaction coefficient" << std::endl;
    std::cout << "  -x|--nx                         Number of points in x direction for SPE" << std::endl;
    std::cout << "  -y|--ny                         Number of points in y direction for SPE" << std::endl;
    std::cout << "  -z|--nz                         Number of points in z direction for SPE" << std::endl;
    std::cout << "  -o|--solver={cheb|pcg}          Outer solver type" << std::endl;
    std::cout << "  -t|--use_tails={yes|no}         Do not use tail removing" << std::endl;
    std::cout << "  -S|--optimize-storage={yes|no}  Do not optimize storage for symmetric matrices" << std::endl;
    std::cout << "  -m|--matrix                     Matrix input file" << std::endl;
    std::cout << "  -a|--ntests                     Number of tests to perform" << std::endl;
    std::cout << "  -p|--prec={uh_cheb|amg|diag}    Preconditioner type" << std::endl;
}

int set_params(int argc, char * argv[], Config& cfg) {
    cfg.niters.push_back(2);
    cfg.sigmas.push_back(3);

    cfg.c                = 1.;
    cfg.nx               = 60;
    cfg.ny               = 220;
    cfg.nz               = 85;

    cfg.ntests           = 1;
    cfg.matrix           = std::string("matrix.dat");
    cfg.use_tails        = 1;
    cfg.optimize_storage = 1;
    cfg.solver           = PCG_SOLVER;
    cfg.prec             = UH_CHEB_PREC;

    static struct option long_options[] = {
	{"sigmas",		required_argument,  NULL, 's'},
	{"niters",		required_argument,  NULL, 'b'},
	{"help",		no_argument,	    NULL, 'h'},
	{"c",			required_argument,  NULL, 'c'},
	{"nx",			required_argument,  NULL, 'x'},
	{"ny",			required_argument,  NULL, 'y'},
	{"nz",			required_argument,  NULL, 'z'},
	{"solver",		required_argument,  NULL, 'o'},
	{"use-tails",		required_argument,  NULL, 't'},
	{"optimize-storage",	required_argument,  NULL, 'S'},
	{"matrix",		required_argument,  NULL, 'm'},
	{"ntests",		required_argument,  NULL, 'a'},
	{"prec",		required_argument,  NULL, 'p'},
	{0, 0, 0, 0}
    };
    while (1) {
	int option_index = 0;
	int ch = getopt_long(argc, argv, "hb:s:S:t:c:x:y:z:to:m:a:p:", long_options, &option_index);

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
	    case 'S': if (!strcmp(optarg, "yes"))
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
		      else
			  THROW_EXCEPTION("Unknown solver type \"" << optarg << "\"");
		      break;
	    case 'm': cfg.matrix = std::string(optarg); break;
	    case 'a': cfg.ntests = uint(atoi(optarg)); break;
	    case 'p': if (!strcmp(optarg, "uh_cheb"))
			  cfg.prec = UH_CHEB_PREC;
		      else if (!strcmp(optarg, "amg"))
			  cfg.prec = AMG_PREC;
		      else if (!strcmp(optarg, "diag"))
			  cfg.prec = DIAG_PREC;
		      else
			  THROW_EXCEPTION("Unknown solver type \"" << optarg << "\"");
		      break;
	    case '?': 
	    default:
		      abort();
	}
    }

    if (cfg.niters.size() == 1 && cfg.sigmas.size() == 1 &&
	((cfg.niters[0] == 2 && cfg.sigmas[0] >= 4) ||
	 (cfg.niters[0] == 3 && cfg.sigmas[0] >= 9))) {
	LOG_WARN("Possibly wrong relation between niter and sigma");
    }

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
