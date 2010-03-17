#include "multi_split_prec.h"
#include "multi_split_misc.h"

#include "include/time.h"
#include "include/tools.h"
#include "include/uvector.h"
#include "include/logger.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>

DEFINE_LOGGER("MultiSplitPrec");

typedef class LinkTypeMultiSplit LinkType;

/*
 * Insertion sort of vector sorted with respect to absolute values in vector a
 * NOTE: later if for some matrices we would get many elements in a row we could
 * replace this sort with a faster one (think heapsort)
 * T must support fabs and < operators
 */
static void psort(const double *a, uint n, uvector<uint>& sorted) {
    sorted[0] = 0;
    double v;
    int j;
    for (uint i = 1; i < n; i++) {
	v = fabs(a[i]);
	for (j = i-1; j >= 0 && fabs(a[sorted[j]]) > v; j--)
	    sorted[j+1] = sorted[j];
	sorted[j+1] = i;
    }
}

void MultiSplitPrec::construct_level(uint level, const SkylineMatrix& A) {
    Level& li = levels[level];
    Level& ln = levels[level+1];

    uvector<double>& aux = li.aux;
    uvector<uint>& tr	 = li.tr;
    uvector<uint>& dtr	 = li.dtr;

    li.N   = A.size();
    li.nnz = A.ja.size();

    /* Number of nodes in the subdomain */
    uint N = li.N;

    uvector<int> nlinks_out(N);
    uvector<int> nlinks_in(N, 0);
    LinkType ltype(A);

    /* Marking stage */
    uint MAX_NUM = 1000; /* Maximum number of elements in a row */
    uvector<uint> sorted(MAX_NUM);
    const double* adata = &(A.a[0]);
    for (uint i = 0; i < N; i++) {
	uint rstart = A.ia[i];		    /* Row start */
	uint rend   = A.ia[i+1];	    /* Row end */
	uint nrz    = rend - rstart - 1;    /* Number of outgoing links */

	nlinks_out[i] = nrz;

	ASSERT(nrz <= MAX_NUM, "Number of nonzero elements in a row = " << nrz);

	for (uint _j = rstart+1; _j < rend; _j++) {
	    uint j = A.ja[_j];
	    nlinks_in[j]++;
	}

	/* Sort off-diagonal elements wrt their abs values */
	psort(adata + rstart+1, nrz, sorted);

	double s = li.q/(1 - li.q) * aux[i];
	for (uint k = 0; k < nrz; k++) {
	    uint _j = rstart+1 + sorted[k];

	    double aij = -A.a[_j];
	    if (aij <= s) {
		uint j = A.ja[_j];

		ltype.mark(i,j); /* mark outgoing link as removable */
		nlinks_out[i]--;
		nlinks_in[j]--;

		s -= aij;
	    } else {
		/* We exhausted available value of c. No other adjoint links can be removed */
		break;
	    }
	}
    }

#if 0
    {
	uint out = 0, in = 0;
	for (uint i = 0; i < N; i++) {
	    out += nlinks_out[i];
	    in  += nlinks_in[i];
	}
	if (out != in)
	    THROW_EXCEPTION("out = " << out << ", in = " << in);

    }
#endif

    /* Compute n & nnz */
    uint n = 0, nnz = 0;
    for (uint i = 0; i < N; i++)
	if (nlinks_out[i] > 0 || nlinks_in[i] > 0) {
	    n++;
	    nnz += nlinks_out[i];
	}
    nnz += n;

    /*
     * Construct next level local matrix
     * For simplicity during the construction we use indices from this level and not the next one yet
     */
    SkylineMatrix& nA = ln.A;

    nA.ia.resize(n+1);
    nA.ja.resize(nnz);
    nA.a.resize(nnz);
    tr.resize(n);
    dtr.reserve(N-n);
    ln.aux.resize(n);

    /*
     * Reverse translation vector: indices of level l -> indices of level l+1
     * NOTE: it is a local vector
     */
    uvector<int> lrevtr(N);

    nA.ia[0] = 0;
    uint ind = 0, iaind = 0;
    for (uint i = 0; i < N; i++) {
	if (nlinks_out[i] > 0 || nlinks_in[i] > 0) {
	    /* The node goes to the next level */

	    /* Calculate the position of the diagonal element */
	    uint dind = ind;

	    double c1 = 1.0/(1 - li.q) * aux[i];

	    nA.ja[dind] = i;
	    nA.a[dind]  = c1;
	    ind++;

	    for (uint _j = A.ia[i]+1; _j < A.ia[i+1]; _j++) {
		uint j = A.ja[_j];
		if (!ltype.is_removed(i,j)) {
		    /* Link goes to coarse level */
		    double v = A.a[_j];

		    nA.ja[ind]  =  j;
		    nA.a[ind]   =  v;
		    nA.a[dind] += -v;

		    ind++;
		}
	    }

	    nA.ia[iaind+1] = ind;
	    tr[iaind] = i;
	    lrevtr[i] = iaind;

	    /* Calculate next level aux array */
	    ln.aux[iaind] = c1;

	    iaind++;

	} else {
	    lrevtr[i] = -1;
	    dtr.push_back(i);
	}
    }
    nA.nrow = nA.ncol = n;

    /* Change matrix to use next level indices */
    for (uint k = 0; k < nA.ja.size(); k++) {
	ASSERT(lrevtr[nA.ja[k]] != -1, "Trying to invert wrong index: k = " << k << ", nA.ja[k] = " << nA.ja[k]);
	nA.ja[k] = lrevtr[nA.ja[k]];
    }

    if (n) {
	/* Allocate space for Chebyshev vectors */
	li.x1.resize(n);
	li.f1.resize(n);
	li.u0.resize(n);
	li.r.resize(n);

	ASSERT(level+1 < nlevels, "Too many levels: " << level+1);
	construct_level(level+1, nA);
    } else {
	LOG_DEBUG("Reducing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level + 1;
    }
}

MultiSplitPrec::MultiSplitPrec(const SkylineMatrix& A, const Config& cfg) : level0_A(A) {
    nlevels = 30;
    levels.resize(nlevels);
    Level& li = levels[0];

    /* Fill in q parameter for each level */
    for (uint l = 0; l < nlevels; l++)
	levels[l].q = (l >= cfg.sigmas.size() ? cfg.sigmas.back() : cfg.sigmas[l]);

    /* Set number of iterations per level */
    for (uint l = 0; l < nlevels; l++)
	levels[l].niter = (l >= cfg.niters.size() ? cfg.niters.back() : cfg.niters[l]);

    uint N = A.size();
    li.aux.resize(N);
    for (uint i = 0; i < N; i++) {
	li.aux[i] = 0.0;
	for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
	    li.aux[i] += A.a[j];
    }

    construct_level(0, A);

    levels.resize(nlevels);
}
