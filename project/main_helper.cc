#include "main.h"
#include <getopt.h>

int set_params(int argc, char * argv[], uint& nlevels, double& c, double& eps, uint& ncheb) {
    nlevels = 12;
    ncheb   = 2;
    eps     = 3.;
    c       = 1.;

    while (1) {
	static struct option long_options[] = {
	    {"nlevels",	required_argument, 0, 'l'},
	    {"eps",	required_argument, 0, 'e'},
	    {"ncheb",	required_argument, 0, 'b'},
	    {"c",	required_argument, 0, 'c'}
	};
	int option_index = 0;
	int c = getopt_long(argc, argv, "l:e:c:", long_options, &option_index);

	if (c == -1)
	    break;

	switch (c) {
	    case 'l': nlevels = uint(atoi(optarg)); break;
	    case 'e': eps = atof(optarg); break;
	    case 'c': c = atof(optarg); break;
	    case 'b': ncheb = uint(atoi(optarg)); break;
	    case '?': 
		      break;
	    default:
		      return 1;
	}
    }
    if (!nlevels || !ncheb || !eps)
	return 2;

    return 0;
}
