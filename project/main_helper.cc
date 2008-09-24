#include "include/logger.h"
#include "main.h"

#include <getopt.h>
#include <cstdlib>

static void usage() {
    std::cout << "Usage: ./main [-e <eps>] [-c <c>] [-b <ncheb>] [-w <nwells>]\n";
}

int set_params(int argc, char * argv[], double& c, double& eps, uint& ncheb, uint& nwells) {
    ncheb   = 2;
    eps     = 3.;
    c       = 1.;
    nwells  = 0;

    static struct option long_options[] = {
	{"eps",	    required_argument, NULL, 'e'},
	{"ncheb",   required_argument, NULL, 'b'},
	{"c",	    required_argument, NULL, 'c'},
	{"nwells",  required_argument, NULL, 'w'},
	{"help",          no_argument, NULL, 'h'},
	{0, 0, 0, 0}
    };
    while (1) {
	int option_index = 0;
	int ch = getopt_long(argc, argv, "e:b:c:w:h", long_options, &option_index);

	if (ch == -1)
	    break;

	switch (ch) {
	    case 'e': eps = atof(optarg); break;
	    case 'c': c = atof(optarg); break;
	    case 'b': ncheb = uint(atoi(optarg)); break;
	    case 'w': nwells = uint(atoi(optarg)); break;
	    case 'h': usage();
		      exit(0);
	    case '?': 
	    default:
		      abort();
	}
    }

    if (!ncheb || !eps) 
	return 2;
    if (ncheb == 2 && eps >= 4 ||
	ncheb == 3 && eps >= 9) {
	std::cout << "Wrong relation between ncheb and eps" << std::endl;
	return 3;
    }


    return 0;
}
