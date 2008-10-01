#include "include/logger.h"
#include "main.h"

#include <getopt.h>
#include <cstdlib>
#include <algorithm>
#include <iterator>

DEFINE_LOGGER("Main");

static void usage() {
    std::cout << "Usage: ./main [-e <sigma>] [-s \"<sigmas>\"] [-c <c>] [-b <niter>] [-w <nwells>]\n";
}

int set_params(int argc, char * argv[], double& c, double& sigma, std::vector<double>& sigmas, 
	       uint& niter, uint& nwells) {
    niter   = 2;
    sigma   = 3.;
    c       = 1.;
    nwells  = 0;

    static struct option long_options[] = {
	{"eps",     required_argument, NULL, 'e'},
	{"sigmas",  required_argument, NULL, 's'},
	{"niter",   required_argument, NULL, 'b'},
	{"c",	    required_argument, NULL, 'c'},
	{"nwells",  required_argument, NULL, 'w'},
	{"help",          no_argument, NULL, 'h'},
	{0, 0, 0, 0}
    };
    while (1) {
	int option_index = 0;
	int ch = getopt_long(argc, argv, "e:s:b:c:w:h", long_options, &option_index);

	if (ch == -1)
	    break;

	switch (ch) {
	    case 'e': sigma = atof(optarg); break;
	    case 's': { 
			  std::istringstream os(optarg);
			  std::copy(std::istream_iterator<double>(os), std::istream_iterator<double>(), 
				    std::back_inserter(sigmas));
		      }
		      break;
	    case 'c': c = atof(optarg); break;
	    case 'b': niter = uint(atoi(optarg)); break;
	    case 'w': nwells = uint(atoi(optarg)); break;
	    case 'h': usage();
		      exit(0);
	    case '?': 
	    default:
		      abort();
	}
    }

    if (!niter || !sigma) 
	return 2;
    if (niter == 2 && sigma >= 4 ||
	niter == 3 && sigma >= 9) {
	LOG_WARN("Possibly wrong relation between niter and sigma");
    }
    if (sigmas.empty())
	sigmas.push_back(sigma);


    return 0;
}
