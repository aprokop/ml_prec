#include "include/logger.h"
#include "main.h"

#include <getopt.h>
#include <cstdlib>
#include <algorithm>
#include <iterator>

DEFINE_LOGGER("Main");

static void usage() {
    std::cout << "Usage: ./main [-e <sigma>] [-s \"<sigmas>\"] [-c <c>] [-b <niter>] [-t <nthreads>] [-w <nwells>]\n";
    std::cout << "		[-nx <nx>] [-ny <ny>] [-nz <nz>]\n";
}

int set_params(int argc, char * argv[], Config& cfg) {
    cfg.niter    = 2;
    cfg.sigma    = 3.;
    cfg.c        = 1.;
    cfg.nwells   = 0;
    cfg.nthreads = 1;

    cfg.nx = 60; 
    cfg.ny = 220; 
    cfg.nz = 85;

    static struct option long_options[] = {
	{"eps",     required_argument, NULL, 'e'},
	{"sigmas",  required_argument, NULL, 's'},
	{"niter",   required_argument, NULL, 'b'},
	{"c",	    required_argument, NULL, 'c'},
	{"nwells",  required_argument, NULL, 'w'},
	{"help",          no_argument, NULL, 'h'},
	{"nthreads",required_argument, NULL, 't'},
	{"nx",	    required_argument, NULL, 'x'},
	{"ny",	    required_argument, NULL, 'y'},
	{"nz",	    required_argument, NULL, 'z'},
	{0, 0, 0, 0}
    };
    while (1) {
	int option_index = 0;
	int ch = getopt_long(argc, argv, "e:s:b:c:w:ht:", long_options, &option_index);

	if (ch == -1)
	    break;

	switch (ch) {
	    case 'e': cfg.sigma = atof(optarg); break;
	    case 's': { 
#ifndef RELX_PREC
			  LOG_WARN("sigmas are given though type is not RELX_PREC");
#endif
			  std::istringstream os(optarg);
			  std::copy(std::istream_iterator<double>(os), std::istream_iterator<double>(), 
				    std::back_inserter(cfg.sigmas));
		      }
		      break;
	    case 'c': cfg.c = atof(optarg); break;
	    case 'b': cfg.niter = uint(atoi(optarg)); break;
	    case 'w': cfg.nwells = uint(atoi(optarg)); break;
	    case 'h': usage();
		      exit(0);
	    case 't': cfg.nthreads = uint(atoi(optarg)); break;
	    case 'x': cfg.nx = uint(atoi(optarg)); break;
	    case 'y': cfg.ny = uint(atoi(optarg)); break;
	    case 'z': cfg.nz = uint(atoi(optarg)); break;
	    case '?': 
	    default:
		      abort();
	}
    }

    if (!cfg.niter || !cfg.sigma) 
	return 2;
    if ((cfg.niter == 2 && cfg.sigma >= 4) ||
	(cfg.niter == 3 && cfg.sigma >= 9)) {
	LOG_WARN("Possibly wrong relation between niter and sigma");
    }
    if (cfg.sigmas.empty())
	cfg.sigmas.push_back(cfg.sigma);

    if (cfg.nthreads > 16) {
	LOG_WARN("Requested number of threads exceeds max value (16)");
	cfg.nthreads = 16;
    }
    LOG_INFO("Number of threads " << cfg.nthreads);

    return 0;
}

std::ostream& operator<<(std::ostream& os, const Config& cfg) {
    os << "mesh    = " << cfg.nx << " x " << cfg.ny << " x " << cfg.nz << std::endl;
    os << "nthread = " << cfg.nthreads << std::endl;
    os << "c       = " << cfg.c << std::endl;
#ifdef CHEB_PREC
    os << "ncheb   = " << cfg.niter << std::endl;
    os << "sigma   = " << cfg.sigma << std::endl;
#elif defined RELX_PREC
    os << "niter   = " << cfg.niter << std::endl;
    os << "gamma   = " << cfg.sigma << std::endl;
    os << "sigmas  = ";
    for (uint i = 0; i < cfg.sigmas.size(); i++)
	os << cfg.sigmas[i] << " ";
    os << std::endl;
#endif
    if (cfg.nwells)
	os << "nwells  = " << cfg.nwells << std::endl;
    return os;
}
