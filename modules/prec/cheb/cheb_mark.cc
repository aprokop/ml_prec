#include "cheb_prec.h"
#include "cheb_misc.h"

#include "include/logger.h"

DEFINE_LOGGER("Prec");

void Prec::construct_permutation(const SkylineMatrix& A, const LinkTypeCheb& ltype_, uvector<int>& nlinks,
				 uint& Md, uint& M, uvector<uint>& map) const {
    uint N = A.size();

    /*
     * We make a copy 'cause we might modify during the construction
     * Takes about 2-3% of total time
     */
    LinkTypeCheb ltype(ltype_);

    uint pind = 0;

    /* Set the marking of all nodes to default */
    uvector<uint> marked(N, 0);

    /* Group 1: all tails */
    if (use_tails) {
	uint i0 = -1, i1 = -1;

	for (uint i = 0; i < N; i++)
	    if (nlinks[i] == 1) {
		/* Start tail with node i */
		i0 = i;

		do {
		    /* Mark node and set new index */
		    map[pind++] = i0;
		    marked[i0]  = 1;

		    /* Find the remaining link */
		    uint j_;
		    ltype.set_row(i0);
		    for (j_ = A.ia[i0]+1; j_ < A.ia[i0+1]; j_++)
			if (ltype.stat(j_) == PRESENT) {
			    i1 = A.ja[j_];
			    break;
			}
		    ASSERT(j_ != A.ia[i0+1], "Could not find remaining link");

		    /* Mark the link as removed */
		    ltype.remove(i0, i1);

		    nlinks[i0] = -1;
		    nlinks[i1]--;

		    i0 = i1;
		} while (nlinks[i0] == 1);

		if (nlinks[i0] == 0) {
		    /* i0 is the end node in fully tridiagonal matrix */
		    map[pind++] = i0;
		    marked[i0]  = 1;

		    nlinks[i0] = -1;
		}
	    }
    }

#if 0
    /* Group 2 : all nodes with two links
     * NOTE: we don't mark links as removed, as we assume that this is the last step */
    const int max_links = 2;
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0 && nlinks[i] <= max_links) {
	    map[pind++] = i;
	    nlinks[i] = -1;

	    marked[i] = 1;
	}
#endif

    M = pind;

    /* Mark the diagonal nodes last */
    uint eind = N-1;
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0) {
	    if (nlinks[i] == 0)
		map[eind--] = i;
	    else
		map[pind++] = i;
	}

    Md = N - pind;
}
