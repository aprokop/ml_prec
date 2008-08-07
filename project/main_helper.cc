#include "main.h"
#include <getopt.h>

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
	{0, 0, 0, 0}
    };
    while (1) {
	int option_index = 0;
	int ch = getopt_long(argc, argv, "e:b:c:w:", long_options, &option_index);

	if (ch == -1)
	    break;

	switch (ch) {
	    case 'e': eps = atof(optarg); break;
	    case 'c': c = atof(optarg); break;
	    case 'b': ncheb = uint(atoi(optarg)); break;
	    case 'w': nwells = uint(atoi(optarg)); break;
	    case '?': 
		      break;
	    default:
		      return 1;
	}
    }

    if (!ncheb || !eps) 
	return 2;

    return 0;
}
