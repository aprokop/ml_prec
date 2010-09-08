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

static uvector<double> aux;
void MultiSplitPrec::construct_level(uint level, const SkylineMatrix& A) {
    Level& li = levels[level];
    Level& ln = levels[level+1];

    uvector<uint>& tr        = li.tr;
    uvector<uint>& dtr       = li.dtr;
    uvector<double>& dtr_val = li.dtr_val;

    li.N   = A.size();
    li.nnz = A.ja.size();

    /* Number of nodes in the subdomain */
    uint N = li.N;

    uvector<int> nlinks_out(N);
    uvector<int> nlinks_in(N, 0);
    LinkType ltype(A);

    /* Marking stage */
    uint MAX_NUM = 100; /* Maximum number of elements in a row */
    uvector<uint> sorted(MAX_NUM);
    const double* adata = &(A.a[0]);
    for (uint i = 0; i < N; i++) {
	uint rstart = A.ia[i];		    /* Row start */
	uint rend   = A.ia[i+1];	    /* Row end */
	uint nrz    = rend - rstart - 1;    /* Number of outgoing links */

	nlinks_out[i] = nrz;

	ASSERT(nrz <= MAX_NUM, "Number of nonzero elements in a row = " << nrz);

	for (uint j_ = rstart+1; j_ < rend; j_++) {
	    uint j = A.ja[j_];
	    nlinks_in[j]++;
	}

	/* Sort off-diagonal elements wrt their abs values */
	psort(adata + rstart+1, nrz, sorted);

	double s = li.q/(1 - li.q) * aux[i];
	for (uint k = 0; k < nrz; k++) {
	    uint j_ = rstart+1 + sorted[k];

	    double aij = -A.a[j_];
	    if (aij <= s) {
		uint j = A.ja[j_];

		ltype.remove(j_); /* mark outgoing link as removable */
		nlinks_out[i]--;
		nlinks_in[j]--;

		s -= aij;
	    } else {
		/*
		 * We exhausted available value of c.
		 * No other adjoint links can be removed
		 * One could split (if s > 0) but we don't do it now
		 */
		break;
	    }
	}

	/* Update c value */
	aux[i] *= 1./(1-li.q);
    }

    std::vector<Tail>& tails = li.tails;
    if (use_tails) {
	uint i0 = uint(-1), i1 = uint(-1);
	for (uint i = 0; i < N; i++) {
	    if (nlinks_out[i] == 1 && nlinks_in[i] <= 1) {
		Tail tail;
		tail.reserve(5);
		TailNode tn;

		double v = 0;

		i0 = i;
		do {
		    /* Step 1: check that there is only one outgoing link */
		    if (nlinks_out[i0] != 1)
			break;

		    /* Step 2: check that there is either no incoming connections,
		     *	       or one incoming connection with the same point, as
		     *	       outgoing
		     */
		    uint j_;
		    for (j_ = A.ia[i0]+1; j_ < A.ia[i0+1]; j_++)
			if (ltype.stat(j_) == PRESENT)
			    break;
		    if (j_ == A.ia[i0+1])
			THROW_EXCEPTION("Row " << i0 << ": remaining link was not found");

		    i1 = A.ja[j_];

		    bool aji_present = (ltype.stat(i1,i0) == PRESENT);
		    if (!( nlinks_in[i0] == 0 ||
			  (nlinks_in[i0] == 1 && aji_present))) {
			/*
			 * Node i0 has more than one incoming connection, or
			 * it has one incoming connection but not from i1
			 */
			break;
		    }

		    ltype.remove(j_);

		    tn.index = i0;
		    tn.a3 = -v; /* Old value */

		    double aij = A.a[j_];   /* < 0 */

		    tn.a2  = 1/(aux[i0] - aij);
		    tn.a1  = -aij*tn.a2;
		    tn.a3 *= tn.a2;

		    double aji = 0;
		    if (aji_present) {
			aji = A(i1,i0);
			ltype.remove(i1, i0);

			nlinks_in[i0]--;
			nlinks_out[i1]--;
		    }
		    v = aji;

		    nlinks_out[i0] = -1;
		    nlinks_in[i1]--;

		    aux[i1] += aux[i0]*(-aji)/(aux[i0] - aij);

		    i0 = i1;

		    tail.push_back(tn);
		} while (true);

		if (!tail.size())
		    continue;

		tn.index = i0;
		tn.a1 = 0.0;

		if (nlinks_out[i0] == 0 && nlinks_in[i0] == 0) {
		    /* i0 is the end node in fully tridiagonal matrix */
		    tn.a2 = 1/aux[i0];
		    nlinks_out[i0] = -1;
		    tail.end_type = 'f';
		} else {
		    /* i0 has incoming/outgoing links */
		    tn.a2 = 1.0;
		    tail.end_type = 'l';
		}
		tn.a3 = -v*tn.a2;

		tail.push_back(tn);

		tails.push_back(tail);
	    }
	}
    }

    /* Compute n & nnz */
    uint n = 0, nnz = 0;
    for (uint i = 0; i < N; i++)
	if (nlinks_out[i] > 0 || nlinks_in[i] > 0) {
	    n++;
	    nnz += nlinks_out[i]+1;
	}

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
    dtr_val.reserve(N-n);

    /* Reverse translation vector: indices of level l -> indices of level l+1 */
    uvector<int> lrevtr(N);

    nA.ia[0] = 0;
    uint ind = 0, iaind = 0;
    for (uint i = 0; i < N; i++) {
	if (nlinks_out[i] > 0 || nlinks_in[i] > 0) {
	    /* The node goes to the next level */

	    /* Save the position of the diagonal element */
	    uint dind = ind;

	    nA.ja[dind] = i;
	    nA.a[dind]  = aux[i];
	    ind++;

	    for (uint j_ = A.ia[i]+1; j_ < A.ia[i+1]; j_++) {
		uint j = A.ja[j_];
		if (ltype.stat(j_) == PRESENT) {
		    /* Link goes to coarse level */
		    double v = A.a[j_];

		    nA.ja[ind]  =  j;
		    nA.a[ind]   =  v;
		    nA.a[dind] += -v;

		    ind++;
		}
	    }

	    nA.ia[iaind+1] = ind;
	    tr[iaind] = i;
	    lrevtr[i] = iaind;

	    /* Calculate new c value */
	    /* NOTE:
	     * At this point aux is part old/part new:
	     * 0, ..., iaind  : new aux
	     * iaind+1, ..., N: old aux
	     */
	    aux[iaind] = aux[i];

	    iaind++;

	} else {
	    lrevtr[i] = -1;
	    if (nlinks_out[i] != -1) {
		/* Link does not belong to a tail */
		dtr.push_back(i);
		dtr_val.push_back(1.0/aux[i]);
	    }
	}
    }
    nA.nrow = nA.ncol = n;
    aux.resize(n);

    /* Change matrix to use next level indices */
    for (uint k = 0; k < nA.ja.size(); k++) {
	ASSERT(lrevtr[nA.ja[k]] != -1, "Trying to invert wrong index: k = " << k << ", nA.ja[k] = " << nA.ja[k]);
	nA.ja[k] = lrevtr[nA.ja[k]];
    }

    if (use_tails) {
	/* Check whether the ends are really local (or belong to a T intersection) */
	for (uint i = 0; i < tails.size(); i++) {
	    Tail& tail = tails[i];
	    if (tail.end_type == 'l' && lrevtr[tail.back().index] == -1)
		tail.end_type = 't';
	}
    }

    li.u0.resize(N);
    li.r.resize(N);
    if (n) {
	/* Allocate space for Chebyshev vectors */
	li.r1.resize(n);
	li.u1.resize(n);

	if (level+1 >= nlevels)
	    THROW_EXCEPTION("Too many levels: " << level+1);

	construct_level(level+1, nA);
    } else {
	LOG_DEBUG("Reducing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level + 1;
    }
}

MultiSplitPrec::MultiSplitPrec(const SkylineMatrix& A, const Config& cfg) : level0_A(A) {
    use_tails = cfg.use_tails;

    nlevels = 30;
    levels.resize(nlevels);

    /* Fill in q parameter for each level */
    for (uint l = 0; l < nlevels; l++)
	levels[l].q = (l >= cfg.sigmas.size() ? cfg.sigmas.back() : cfg.sigmas[l]);

    /* Set number of iterations per level */
    levels[0].niter = 1;
    for (uint l = 1; l < nlevels; l++)
	levels[l].niter = (l >= cfg.niters.size() ? cfg.niters.back() : cfg.niters[l]);

    uint N = A.size();
    aux.resize(N);
    for (uint i = 0; i < N; i++) {
	aux[i] = 0.0;
	for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
	    aux[i] += A.a[j];
    }

    construct_level(0, A);

    levels.resize(nlevels);
}
