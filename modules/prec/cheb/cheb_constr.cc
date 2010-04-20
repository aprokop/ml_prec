#include "cheb_prec.h"
#include "cheb_misc.h"

#include "include/time.h"
#include "include/tools.h"
#include "include/uvector.h"
#include "include/svector.h"
#include "include/logger.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <set>

DEFINE_LOGGER("Prec");

typedef class LinkTypeCheb LinkType;

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

    li.N   = A.size();
    li.nnz = A.ja.size();

    /* Number of nodes in the subdomain */
    uint N = li.N;

    uvector<int> nlinks(N, 0);
    LinkType ltype(A);

    uvector<double>& aux = li.aux;
    aux.resize(N);

    /* Drop links */
    log_state("M");
    uint MAX_NUM = 1000; /* Maximum number of elements in a row */
    uvector<uint> sorted(MAX_NUM);
    const double* adata = A.a.data();
    for (uint i = 0; i < N; i++) {
	uint rstart = A.ia[i];		    /* Row start */
	uint rend   = A.ia[i+1];	    /* Row end */
	uint nrz    = rend - rstart - 1;    /* Number of links */

	ASSERT(nrz <= MAX_NUM, "Number of nonzero elements in a row = " << nrz);

	nlinks[i] += nrz;

	/* TODO: do not compute aux[i] every level; translate from upper to lower */
	aux[i] = 0.0;
	for (uint j_ = rstart; j_ < rend; j_++)
	    aux[i] += A.a[j_];

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
		if (ltype.mark(i,j) == REMOVED) {
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
    log_state("d");

    uint& M             = li.M;
    uint& Md		= li.Md;
    uvector<uint>& map  = li.map;
    uvector<uint>& rmap = li.rmap;

    map.resize(N);
    log_state("P");
    construct_permutation(A, ltype, nlinks, Md, M, map);
    log_state("d");

    /* Construct reverse map */
    rmap.resize(N);
    for (uint i = 0; i < N; i++)
	rmap[map[i]] = i;

    SkylineMatrix& nA = ln.A;
    SkylineMatrix&  U = li.U;
    CSRMatrix&      L = li.L;

    /* Saad. Iterative methods for sparse linear systems. Pages 310-312 */
    uvector<int> jr(N, -1);
    // typedef std::set<uint> container;
    typedef svector<uint> container;
    container jw;
    uvector<double> w;
    uint max_num;

    uint n = N-M-Md;
    L.ia.push_back(0);     L.nrow = N-Md;    L.ncol = N-Md;
    U.ia.push_back(0);     U.nrow = M;       U.ncol = N-Md;
    nA.ia.push_back(0);   nA.nrow = n;      nA.ncol = n;

    /* Reserve some space */
    // L.ja.reserve(2*(N-Md));       L.a.reserve(2*(N-Md));
    // U.ja.reserve(3*M);		  U.a.reserve(3*M);
    // nA.ja.reserve(5*(N-M-Md));    nA.a.reserve(5*(N-M-Md));

    /* TODO: deal with M = 0 */
    /* i corresponds to a permuted index */
    log_state("LU");
    for (uint i = 0; i < N-Md; i++) {
	/* Step 0: clear tmp values */
	unsigned jwn = jw.size();
	for (uint k = 0; k < jwn; k++)
	    jr[jw[k]] = -1;
	jw.clear();
	w.clear();

	uint arow = map[i];	/* Index of a row in A */

	/* Step 1: create buffer with permuted row of A */
	/* Add diagonal element to buffer */
	jr[i] = 0;
	jw.insert(i);
	w.push_back(aux[map[i]]);    /* w[0] is the value of the diagonal element */

	/* Add off-diagonal elements to buffer */
	max_num = 1;
	for (uint j_ = A.ia[arow]+1; j_ < A.ia[arow+1]; j_++) {
	    uint j = A.ja[j_];

	    if (ltype.stat(arow,j) == PRESENT) {
		uint new_j = rmap[j];   /* permuted index */

		jr[new_j] = max_num++;
		jw.insert(new_j);

		double z = A.a[j_] / li.beta;   /* scale the element */

		w.push_back(z);
		w[0] -= z;       /* update diagonal */
	    }
	}

	/* Step 2: perform sparse gaussian elimination */
	uint m = std::min(i, M);

	uint k = *(jw.begin());
	while (k < m) {
	    double lik = w[jr[k]] / U(k,k);

	    /* Update L */
	    L.ja.push_back(k);
	    L.a.push_back(lik);

	    /* Update buffer using k-th row of U */
	    for (uint j_ = U.ia[k]+1; j_ < U.ia[k+1]; j_++) {
		uint j = U.ja[j_];      /* j is a permuted index */

		if (jr[j] != -1) {
		    /* Element already exists in the buffer => update */
		    w[jr[j]] -= lik*U.a[j_];
		} else {
		    /* Element does not exist in the buffer => create */
		    jr[j] = max_num++;
		    jw.insert(j);
		    w.push_back(-lik*U.a[j_]);
		}
	    }

	    k = *(jw.upper_bound(k));
	}
	L.ia.push_back(L.ja.size());

	/* Step 3: move element from buffer to corresponding rows of U/A_{level+1} */
	if (i < M) {
	    /* Update U */
	    for (container::const_iterator it = jw.lower_bound(i); it != jw.end(); it++) {
		U.ja.push_back(*it);
		U.a.push_back(w[jr[*it]]);
	    }

	    U.ia.push_back(U.ja.size());
	} else {
	    /* Update A_{level+1} */
	    uint adind = nA.ja.size();
	    nA.ja.push_back(i-M);
	    nA.a.push_back(w[0]);

	    for (container::const_iterator it = jw.lower_bound(M); it != jw.end(); it++)
		if (*it != i) {
		    nA.ja.push_back(*it-M);
		    nA.a.push_back(w[jr[*it]]);
		}

	    nA.ia.push_back(nA.ja.size());
	}
    }
    log_state("d");

    /* Process diagonal block */
    log_state("D");
    uvector<double>& dval = li.dval;
    dval.resize(Md);
    for (uint i = 0; i < Md; i++)
	dval[i] = 1./aux[map[i+(N-Md)]];
    log_state("d");

    li.w.resize(N);

    if (n) {
	/* Allocate space for Chebyshev vectors */
	li.tmp.resize(n);
	li.x2.resize(n);
	li.u0.resize(n);
	li.u1.resize(n);

	if (level+1 >= nlevels)
	    THROW_EXCEPTION("Too many levels: " << level+1);

	construct_level(level+1, nA);
    } else {
	LOG_DEBUG("Reducing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level + 1;
    }
}

Prec::Prec(const SkylineMatrix& A, const Config& cfg) : level0_A(A) {
    /* Initialize ViTE */
    foss = new std::ostringstream;
    (*foss) << std::fixed << std::setprecision(8);
    vite_start = pclock();

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
	/* Symmetrize the matrix */
	/* NOTE: level0_A reference WILL BE incorrect, as it references the unsymmetric matrix */
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

Prec::~Prec() {
    dump_vite_trace();
}

/* Optimize level matrices for symmetricity */
void Prec::optimize_storage() {
    for (uint l = 1; l < nlevels; l++)
	levels[l].A.optimize_storage();
}
