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

    /* Initialize LinkType; the status of each link is PRESENT */
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
        ltype.set_row(i);
        for (uint _j = 0; _j < nrz; _j++) {
            double v = A.a[rstart+1 + sorted[_j]];

            if (to_remove(aux[i], v, li.beta, s)) {
                /* It is possible to remove the link from i-th end */
                uint j_ = rstart+1 + sorted[_j];
                uint j = A.ja[j_];

                /* Mark the link for removal from i-th end; the new link status is returned */
                if (ltype.mark(j_) == REMOVED) {
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
    /*
     * Construct node perumutation.
     * All nodes are divided into three groups:
     *	0, ..., M-1   : Nodes which are connected to some nodes and which will be excluded
     *  M, ..., N-Md  : Nodes which go to the next level
     *  Md, ..., N    : Nodes which have no connections to other nodes (diagonal submatrix). Excluded
     */
    construct_permutation(A, ltype, nlinks, Md, M, map);
    log_state("d");

    /* Construct reverse map */
    rmap.resize(N);
    for (uint i = 0; i < N; i++)
        rmap[map[i]] = i;

    SkylineMatrix& nA = ln.A;
    SkylineMatrix&  U = li.U;
    CSRMatrix&      L = li.L;

    log_state("LU");
    construct_sparse_lu(A, map, rmap, Md, M, ltype, li.beta, aux, nA, U, L);
    log_state("d");

    /* Process diagonal block */
    log_state("D");
    uvector<double>& dval = li.dval;
    dval.resize(Md);
    for (uint i = 0; i < Md; i++)
        /* Instead of keeping diagonal, keep its reciprocal */
        dval[i] = 1./aux[map[i+(N-Md)]];
    log_state("d");

    li.w.resize(N);

    uint n = nA.size();
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

    /* Initialize stream for dumping norms */
    norm_oss = new std::ostringstream;
    (*norm_oss) << std::scientific;

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
    if (!cfg.nonsym_matrix) {
        /* Given matrix is described as symmetric. Do nothing */
        Asym = const_cast<SkylineMatrix*>(&A);
    } else {
        /*
         * "Symmetrize" the matrix
         * That means that we set Asym(i,j) = max(A(i,j), A(j,i))
         * In other words we increase some offdiagonal elements to achive symmtery. Theory
         * says that such action results in the M-matrix if the original matrix was an M-matrix
         * Our hope is that the preconditioner for the symmetrized matrix would be a good
         * preconditioner for the original matrix
         */
        /* NOTE: level0_A reference WILL BE incorrect, as it references the nonsymmetric matrix */
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

    /* Construct preconditioner */
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

    if (cfg.nonsym_matrix)
        delete Asym;
}

Prec::~Prec() {
    /* Output ViTE trace */
    dump_vite_trace();
    delete foss;

    dump_norm_trace();
    delete norm_oss;
}

/* Optimize level matrices for symmetricity */
void Prec::optimize_storage() {
    for (uint l = 1; l < nlevels; l++)
        levels[l].A.optimize_storage();
}
