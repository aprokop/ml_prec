#include "solvers.h"

#include <cstdlib>

void generate_x0(Vector& x) {
#if 0
    srandom(time(NULL)); 
#else
    // we don't want it to change from run to run
    srandom(3);
#endif
    for (uint i = 0; i < x.size(); i++) {
	x[i] = 20.*(random() - 0.5*RAND_MAX)/RAND_MAX + 100;
	// x[i] = 1 - (i&1)*2;
	// x[i] = 1;
    }
}

