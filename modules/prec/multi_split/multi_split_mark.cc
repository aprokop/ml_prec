#include "multi_split_prec.h"
#include "multi_split_misc.h"

#include "include/logger.h"

DEFINE_LOGGER("Prec");

void MultiSplitPrec::order_none    (const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
				    const uvector<int>& nlinks_in, const uvector<int>& nlinks_out,
				    uint& Md, uint& M, uvector<uint>& map) const {
    uint N = A.size();

    uint pind = 0;
    M = 0;

    /* Mark the diagonal nodes last */
    uint eind = N-1;
    for (uint i = 0; i < N; i++)
	if (nlinks_in[i] == 0 && nlinks_out[i] == 0)
	    map[eind--] = i;
	else
	    map[pind++] = i;

    Md = N - pind;
}

void MultiSplitPrec::order_original(const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
				    uvector<int>& nlinks_in, uvector<int>& nlinks_out,
				    uint& Md, uint& M, uvector<uint>& map) const {
    uint N = A.size();

    /*
     * We make a copy 'cause we might modify during the construction
     * Takes about 2-3% of total time
     */
    LinkTypeMultiSplit ltype(ltype_);

    uint pind = 0;

    /* Set the marking of all nodes to default */
    uvector<uint> marked(N, 0);

    /* Group 1: all tails */
    if (use_tails) {
	uint i0 = -1, i1 = -1;

	for (uint i = 0; i < N; i++)
	    if (nlinks_out[i] == 1 && nlinks_in[i] <= 1) {
		/* Start tail with node i */
		i0 = i;

		do {
		    /* Step 1: check that there is only one outgoing link */
		    if (nlinks_out[i0] != 1)
			break;

		    /* Step 2: find the remaining link */
		    uint j_;
		    for (j_ = A.ia[i0]+1; j_ < A.ia[i0+1]; j_++)
			if (ltype.stat(j_) == PRESENT)
			    break;
		    if (j_ == A.ia[i0+1])
			THROW_EXCEPTION("Row " << i0 << ": remaining link was not found");
		    i1 = A.ja[j_];

		    /* Step 3: check that there is either no incoming connections,
		     *	       or one incoming connection with the same point, as
		     *	       outgoing
		     */
		    bool aji_present = (ltype.stat(i1,i0) == PRESENT);
		    if (!( nlinks_in[i0] == 0 ||
			   (nlinks_in[i0] == 1 && aji_present))) {
			/*
			 * Node i0 has more than one incoming connection, or
			 * it has one incoming connection but not from i1
			 */
			break;
		    }

		    /* Everything is fine. Mark the node and set new index */
		    map[pind++] = i0;
		    marked[i0]  = 1;

		    ltype.remove(j_);

		    if (aji_present) {
			ltype.remove(i1, i0);

			nlinks_in[i0]--;
			nlinks_out[i1]--;
		    }

		    nlinks_out[i0] = -1;
		    nlinks_in[i1]--;

		    i0 = i1;
		} while (true);

		if (nlinks_out[i0] == 0 && nlinks_in[i0] == 0) {
		    /* i0 is the end node in fully tridiagonal matrix */
		    map[pind++] = i0;
		    marked[i0]  = 1;

		    nlinks_out[i0] = -1;
		}
	    }
    }

#if 0
    /* Group 2 : all nodes with very small c */
    uint c_marked = 0;
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0) {
	    double s = 0.0;
	    for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++)
		if (ltype.stat(j) == PRESENT)
		    s -= A.a[j];

	    if (1 - aux[i] / (aux[i] + s) > 0.1) {
		map[pind++] = i;

		nlinks_in[i] = 0;
		nlinks_out[i] = -1;

		marked[i] = 1;
		c_marked++;
	    }
	}
    LOG_INFO("Excluded using small c criteria: " << c_marked);
#endif

#if 0
    /* Group 2 : all nodes with two links
     * NOTE: we don't mark links as removed, as we assume that this is the last step */
    const int max_links = 2;
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0 && (nlinks_in[i] <= max_links || nlinks_out[i] <= max_links)) {
	    map[pind++] = i;
	    nlinks_in[i] = 0;
	    nlinks_out[i] = -1;

	    marked[i] = 1;
	}
#endif

    M = pind;

    /* Mark the diagonal nodes last */
    uint eind = N-1;
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0) {
	    if (nlinks_in[i] == 0 && nlinks_out[i] == 0)
		map[eind--] = i;
	    else
		map[pind++] = i;
	}

    Md = N - pind;
}

void MultiSplitPrec::order_simple_1(const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
				    const uvector<int>& nlinks_in, const uvector<int>& nlinks_out,
				    uint& Md, uint& M, uvector<uint>& map) const {
    uint N = A.size();

    const int MAX_USEFUL_DEGREE = 6;
    std::vector<uvector<uint> > degrees(MAX_USEFUL_DEGREE+2); // degrees[MAX_USEFUL_DEGREE+1] contains all with >

    for (uint i = 0; i < N; i++)
	if (nlinks_out[i] <= MAX_USEFUL_DEGREE)
	    degrees[nlinks_out[i]].push_back(i);
	else
	    degrees.back().push_back(i);

    const uvector<uint>& list0 = degrees[0];

    uint pind = 0;
    for (uint k = 0; k < list0.size(); k++)
	if (nlinks_in[list0[k]] != 0)
	    map[pind++] = list0[k];
    for (uint i = 1; i < degrees.size(); i++) {
	const uvector<uint>& list = degrees[i];
	for (uint k = 0; k < list.size(); k++)
	    map[pind++] = list[k];

	if (i == degree_max)
	    M = pind;
    }

    /* Mark the diagonal nodes last */
    Md = pind;
    for (uint k = 0; k < list0.size(); k++)
	if (nlinks_in[list0[k]] == 0)
	    map[pind++] = list0[k];
    Md = N - Md;
}

#define BSIZE 5
void MultiSplitPrec::order_block(const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
				 const uvector<int>& nlinks_in, const uvector<int>& nlinks_out,
				 uint& Md, uint& M, uvector<uint>& map) const {
    const uint nx = 60;
    const uint ny = 220;
    const uint nz = 85;

    uint N = A.size();
    if (N != nx*ny*nz)
	THROW_EXCEPTION("Wrong mesh size: " << N << " (expected 1122000)");
    LOG_WARN("Mesh is assumed to be cartesian 3D 60x220x85");

    /* Set the marking of all nodes to default */
    uvector<uint> marked(N, 0);

    uint pind = 0;
#if 0
    for (uint k = 0; k < nz; k++)
	for (uint j = 0; j < ny; j++)
	    for (uint i = 0; i < nx; i++) {
		uint id = k*ny*nx + j*nx + i;
		if ( ((i % BSIZE != 0) && (j % BSIZE != 0) && (k % BSIZE != 0)) &&  /* check that it is an inner node */
		     (nlinks_in[id] || nlinks_out[id]) &&			    /* check that it is not a diagonal node */
		     (marked[id] == 0)) {					    /* check that it is unmarked */

		    map[pind++] = id;
		    marked[id] = 1;
		}
	    }
#else
    for (uint K = 0; K < nz; K += BSIZE)
	for (uint J = 0; J < ny; J += BSIZE)
	    for (uint I = 0; I < nx; I += BSIZE) {
		for (uint k = K+1; k < std::min(K+BSIZE,nz)-1; k++)
		    for (uint j = J+1; j < std::min(J+BSIZE,ny)-1; j++)
			for (uint i = I+1; i < std::min(I+BSIZE,nx)-1; i++) {
			    uint id = k*ny*nx + j*nx + i;
			    if (nlinks_in[id] || nlinks_out[id]) { /* check that it is not a diagonal node */

				map[pind++] = id;
				marked[id] = 1;
			    }
			}
	    }
#endif
    M = pind;

    /* Mark the remaining nodes, diagonal nodes */
    uint eind = N-1;
    for (uint i = 0; i < N; i++)
	if (marked[i] == 0) {
	    if (nlinks_in[i] == 0 && nlinks_out[i] == 0)
		map[eind--] = i;
	    else
		map[pind++] = i;
	}

    Md = N - pind;
}
