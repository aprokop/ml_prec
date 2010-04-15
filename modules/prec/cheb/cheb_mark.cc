#include "cheb_prec.h"
#include "cheb_misc.h"

#include "include/logger.h"

DEFINE_LOGGER("Prec");

void Prec::construct_permutation(const SkylineMatrix& A, LinkTypeCheb ltype, uvector<int>& nlinks,
				 uint& M, uvector<uint>& map, uvector<uint>& rmap) const {
    uint N = A.size();

    uint pind = 0;

    uvector<uint> marked(N, 0);

    /* Group 1: all links with zero links */
    for (uint i = 0; i < N; i++)
	if (nlinks[i] == 0) {
	    rmap[i] = pind;
	    map[pind] = i;
	    pind++;

	    marked[i] = 1;
	}

#if 0
    /* Group 2: all tails */
    if (use_tails) {
	uint i0, i1;

	for (uint i = 0; i < N; i++)
	    if (nlinks[i] == 1) {
		/* Start tail with node i */
		i0 = i;

		do {
		    rmap[i0]   = pind;
		    map[pind]  = i0;
		    pind++;

		    marked[i0] = 1;

		    /* Find the remaining link */
		    uint j_;
		    for (j_ = A.ia[i0]+1; j_ < A.ia[i0+1]; j_++) {
			i1 = A.ja[j_];
			if (ltype.stat(i0, i1) == PRESENT)
			    break;
		    }
		    ASSERT(j_ != A.ia[i0+1], "Could not find remaining link");

		    ltype.remove(i0, i1);

		    nlinks[i0] = -1;
		    nlinks[i1]--;

		    i0 = i1;
		} while (nlinks[i0] == 1);

		if (nlinks[i0] == 0) {
		    /* i0 is the end node in fully tridiagonal matrix */
		    rmap[i0]   = pind;
		    map[pind]  = i0;
		    pind++;

		    marked[i0] = 1;

		    nlinks[i0] = -1;
		}
	    }
    }
#endif

    M = pind;

    /* Mark all others */
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0) {
	    rmap[i]   = pind;
	    map[pind] = i;
	    pind++;

	    marked[i] = 1;
	}
}
