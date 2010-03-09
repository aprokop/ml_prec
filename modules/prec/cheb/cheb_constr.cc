#include "cheb_prec.h"

#include "include/time.h"
#include "include/tools.h"
#include "include/uvector.h"
#include "include/logger.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>

DEFINE_LOGGER("Prec");

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

void Prec::construct_level(uint level, const SkylineMatrix& A) {
    Level& li = levels[level];
    Level& ln = levels[level+1];

    uvector<double>& aux	     = li.aux;
    uvector<uint>& tr		     = li.tr;
    uvector<uint>& dtr		     = li.dtr;

    li.N   = A.size();
    li.nnz = A.ja.size();

    /* Number of nodes in the subdomain */
    uint N = li.N;

    std::vector<int> nlinks(N);
    LinkType ltype(A);

    /* Mark */
    uint MAX_NUM = 1000; /* Maximum number of elements in a row */
    uvector<uint> sorted(MAX_NUM);
    const double* adata = A.a.data();
    for (uint i = 0; i < N; i++) {
	uint rstart = A.ia[i];		    /* Row start */
	uint rend   = A.ia[i+1];	    /* Row end */
	uint nrz    = rend - rstart - 1;    /* Number of links */

	nlinks[i] += nrz;

	ASSERT(nrz <= MAX_NUM, "Number of nonzero elements in a row = " << nrz);

	/* Sort off-diagonal elements wrt their abs values */
	psort(adata + rstart+1, nrz, sorted);

	double s = 0.;
	for (uint _j = 0; _j < nrz; _j++) {
	    double v = A.a[rstart+1 + sorted[_j]];

	    if (to_remove(aux[i], v, li.beta, s)) {
		/* The link is marked as free to remove from i-th end */
		uint j = A.ja[rstart+1 + sorted[_j]];

		/* TODO: for j > i we actually already know the offset in LinkType::a array
		 * so it can be significantly sped up */
		if (ltype.mark(i,j)) {
		    /* This link is marked as removable from both ends so it is removed */
		    nlinks[i]--;
		    nlinks[j]--;
		}
	    } else {
		/* We exhausted available value of c. No other adjoint links can be removed */
		break;
	    }
	}
    }

    std::vector<Tail>& tails = li.tails;
    if (use_tails) {
	uint tnn = 0;
	for (uint i = 0; i < N; i++)
	    if (nlinks[i] == 1)
		tnn++;

	/* Remove tails inside the subdomain */
	tails.resize(tnn+1);
	tnn = 0;

	uint i0 = uint(-1), i1 = uint(-1);
	double v = 0;
	for (uint i = 0; i < N; i++)
	    if (nlinks[i] == 1) {
		Tail& tail = tails[tnn];
		tail.reserve(5);
		TailNode tn;

		i0 = i;
		do {
		    tn.index = i0;
		    tn.a3 = -v; /* old v */

		    /* Find the remaining link */
		    uint _j;
		    if (!ltype.remove_remaining_link(i0, _j)) {
			/* We haven't found the links in the matrix => it is a boundary link. Skip it */
			THROW_EXCEPTION("Huh??");
			break;
		    }

		    i1 = A.ja[_j];
		    /* The link actually belongs to the next level so it is scaled */
		    v = A.a[_j] / li.beta;

		    tn.a2  = 1/(aux[i0] - v);
		    tn.a1  = -v*tn.a2;
		    tn.a3 *= tn.a2;

		    nlinks[i0] = -1;
		    nlinks[i1]--;

		    aux[i1] += aux[i0]*(-v) / (aux[i0] - v);

		    i0 = i1;

		    tail.push_back(tn);
		} while (nlinks[i0] == 1);

		if (!tail.size())
		    continue;

		tn.index = i0;

		if (nlinks[i0] == 0) {
		    /* i0 is the end node in fully tridiagonal matrix */
		    tn.a2 = 1/aux[i0];
		    nlinks[i0] = -1;
		    tail.end_type = 'f';
		} else {
		    tn.a2 = 1.;
		    tail.end_type = 'l';
		}
		tn.a3 = -v*tn.a2;

		tail.push_back(tn);

		tnn++;
	    }
	tails.resize(tnn);
    }

    /* Compute n & nnz. Small problem: nnz might include bnd links */
    uint n = 0, nnz = 0;
    for (uint i = 0; i < N; i++)
	if (nlinks[i] > 0) {
	    n++;
	    nnz += nlinks[i];
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
	if (nlinks[i] > 0) {
	    /* The node goes to the next level */

	    /* Calculate the position of the diagonal element */
	    uint dind = ind;

	    nA.ja[ind] = i;
	    nA.a[ind] = aux[i];
	    ind++;

	    for (uint _j = A.ia[i]+1; _j < A.ia[i+1]; _j++) {
		uint j = A.ja[_j];
		if (!ltype.is_removed(i,j)) {
		    /* Link goes to coarse level */
		    /* Scale the link */
		    double v = A.a[_j] / li.beta;

		    nA.ja[ind]  = j;
		    nA.a[ind]   = v;
		    nA.a[dind] += -v;

		    ind++;
		}
	    }

	    nA.ia[iaind+1] = ind;
	    tr[iaind] = i;
	    lrevtr[i] = iaind;

	    /* Calculate next level aux array */
	    ln.aux[iaind] = aux[i];

	    iaind++;

	} else {
	    lrevtr[i] = -1;
	    if (nlinks[i] != -1) {
		/*
		 * All links for this point are gone and the point does not belong to a tail.
		 * Add it to the diagonal vector
		 */
		dtr.push_back(i);
	    }
	}
    }
    nA.nrow = nA.ncol = n;

    /* Change matrix to use next level indices */
    for (uint k = 0; k < nA.ja.size(); k++) {
	ASSERT(lrevtr[nA.ja[k]] != -1, "Trying to invert wrong index: k = " << k << ", nA.ja[k] = " << nA.ja[k]);
	nA.ja[k] = lrevtr[nA.ja[k]];
    }

    if (use_tails) {
	/* Check whether the ends are really local (or belong to a T intersection */
	for (uint i = 0; i < tails.size(); i++) {
	    Tail& tail = tails[i];
	    if (tail.end_type == 'l' && lrevtr[tail.back().index] == -1)
		tail.end_type = 't';
	}
    }

    if (n) {
	/* Allocate space for Chebyshev vectors */
	li.x1.resize(n);
	li.f1.resize(n);
	li.u0.resize(n);
	li.u1.resize(n);
	li.tmp.resize(n);

	ASSERT(level+1 < nlevels, "Too many levels: " << level+1);
	construct_level(level+1, nA);
    } else {
	LOG_DEBUG("Reducing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level + 1;
    }
}

Prec::Prec(const SkylineMatrix& A, const Config& cfg) : level0_A(A) {
    use_tails = cfg.use_tails;

    nlevels = 30;
    levels.resize(nlevels);
    Level& li = levels[0];

    /* Fill in alpha, beta parameters for each level */
    for (uint l = 0; l < nlevels && l < cfg.sigmas.size(); l++) {
	levels[l].alpha = 1.0;
	levels[l].beta  = cfg.sigmas[l];
    }
    for (uint l = cfg.sigmas.size(); l < nlevels; l++) {
	levels[l].alpha = 1.0;
	levels[l].beta = cfg.sigmas.back();
    }

    SkylineMatrix* Asym;
    if (!cfg.unsym_matrix)
	Asym = const_cast<SkylineMatrix*>(&A);
    else {
	/* NOTE: level0_A reference IS incorrect, as it references
	 * the unsymmetric matrix */
	Asym = new SkylineMatrix(A);

	uint N = Asym->size();
	const uvector<uint>& ia = Asym->get_ia();
	const uvector<uint>& ja = Asym->get_ja();
	uvector<double>&      a = Asym->a;

	/* Create a symmetric variant of the matrix */
	for (uint i = 0; i < N; i++)
	    for (uint j = ia[i]+1; j < ia[i+1]; j++) {
		double d = (*Asym)(ja[j],i); /* <= 0 */
		if (d > a[j])
		    a[j] = d;
	    }
    }

    uint N = Asym->size();
    li.aux.resize(N);
    for (uint i = 0; i < N; i++) {
	li.aux[i] = 0;
	for (uint j = Asym->ia[i]; j < Asym->ia[i+1]; j++)
	    li.aux[i] += Asym->a[j];
    }
    construct_level(0, *Asym);

    /* Set number of Chebyshev iterations per level */
    for (uint l = 0; l < nlevels; l++)
	levels[l].ncheb = (l >= cfg.niters.size() ? cfg.niters.back() : cfg.niters[l]);

    /* Calculate constants of spectral equivalence */
    levels[nlevels-1].lmin = levels[nlevels-1].alpha;
    levels[nlevels-1].lmax = levels[nlevels-1].beta;
    for (int l = nlevels-2; l >= 0; l--) {
	Level& li = levels[l];
	Level& ln = levels[l+1];

	double cs = cheb((ln.lmax + ln.lmin)/(ln.lmax - ln.lmin), li.ncheb);
	li.lmin = li.alpha * (1 - 1/cs);
	li.lmax = li.beta  * (1 + 1/cs);
    }

    levels.resize(nlevels);
    if (cfg.optimize_storage) {
	TIME_INIT();
	TIME_START();
	this->optimize_storage();
	LOG_DEBUG(TIME_INFO("Storage optimization"));
    }

    if (cfg.unsym_matrix)
	delete Asym;
}

/* Optimize level matrices for symmetricity */
void Prec::optimize_storage() {
    for (uint l = 1; l < nlevels; l++)
	levels[l].A.optimize_storage();
}
