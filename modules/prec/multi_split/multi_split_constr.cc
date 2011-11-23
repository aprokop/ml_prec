#include "multi_split_prec.h"
#include "multi_split_misc.h"
#include "modules/config_file/config_file.h"

#include "include/time.h"
#include "include/tools.h"
#include "include/uvector.h"
#include "include/logger.h"

/* Needed for DirectSolver_free */
#include "modules/solvers/solvers.h"

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

    li.N   = A.size();
    li.nnz = A.ja.size();

#ifndef HAVE_UMFPACK
    if (level == nlevels-1)
	THROW_EXCEPTION("Too many levels: " << nlevels);
#endif

    if (li.N <= coarse_n || level == nlevels-1) {
	nlevels  = level + 1;
	coarse_n = li.N;
	return;
    }

    /* Number of nodes in the subdomain */
    uint N = li.N;

    uvector<int> nlinks_out(N);
    uvector<int> nlinks_in(N, 0);
    LinkType ltype(A);

    /* Marking stage */
    uint MAX_NUM = 100; /* Maximum number of elements in a row */
    uvector<uint> sorted(MAX_NUM);
    aux.resize(N);
    const double* adata = &(A.a[0]);
    for (uint i = 0; i < N; i++) {
	uint rstart = A.ia[i];		    /* Row start */
	uint rend   = A.ia[i+1];	    /* Row end */
	uint nrz    = rend - rstart - 1;    /* Number of outgoing links */

	nlinks_out[i] = nrz;

	ASSERT(nrz <= MAX_NUM, "Number of nonzero elements in a row = " << nrz);

	aux[i] = A.a[rstart];
	for (uint j_ = rstart+1; j_ < rend; j_++) {
	    uint j = A.ja[j_];
	    nlinks_in[j]++;

	    aux[i] += A.a[j_];
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

    uint& M             = li.M;
    uint& Md		= li.Md;
    uvector<uint>& map  = li.map;
    uvector<uint>& rmap = li.rmap;

    map.resize(N);
    /*
     * Construct node perumutation.
     * All nodes are divided into three groups:
     *	0, ..., M-1   : Nodes which are connected to some nodes and which will be excluded
     *  M, ..., N-Md  : Nodes which go to the next level
     *  Md, ..., N    : Nodes which have no connections to other nodes (diagonal submatrix). Excluded
     */
    switch (order) {
	case ORDER_NONE    : order_none    (A, ltype, aux, nlinks_in, nlinks_out, Md, M, map); break;
	case ORDER_ORIGINAL: order_original(A, ltype, aux, nlinks_in, nlinks_out, Md, M, map); break;
	case ORDER_BLOCK   : order_block   (A, ltype, aux, nlinks_in, nlinks_out, Md, M, map); break;
	case ORDER_SIMPLE_1: order_simple_1(A, ltype, aux, nlinks_in, nlinks_out, Md, M, map); break;
    }

    /* Construct reverse map */
    rmap.resize(N);
    for (uint i = 0; i < N; i++)
	rmap[map[i]] = i;

    SkylineMatrix& nA = ln.A;
    SkylineMatrix&  U = li.U;
    CSRMatrix&      L = li.L;

    construct_sparse_lu(A, map, rmap, Md, M, ltype, 1.0, aux, nA, U, L);

    /* Process diagonal block */
    uvector<double>& dval = li.dval;
    dval.resize(Md);
    for (uint i = 0; i < Md; i++)
	/* Instead of keeping diagonal, keep its reciprocal */
	dval[i] = 1./aux[map[i+(N-Md)]];

    li.u0.resize(N);
    li.r.resize(N);
    li.w.resize(N);

    const uint n = nA.size();
    if (n) {
	/* Allocate space for Chebyshev vectors */
	li.x2.resize(n);
	li.F.resize(n);

	construct_level(level+1, nA);
    } else {
	LOG_DEBUG("Reducing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level + 1;
    }
}

MultiSplitPrec::MultiSplitPrec(const SkylineMatrix& A, const Config& cfg) : level0_A(A) {
    /* Read parameters from the config file */
    std::string inner_iter_type, lu_type, elim_order_type, elim_stop_type;
    std::string qs_str, niters_str, epss_str;

    std::vector<double> qs;
    std::vector<uint>	niters;
    std::vector<double> epss;

    if (!cfg.prec_conf_file.empty()) {
	LOG_INFO("Prec configuration file provided, ignoring command line arguments");
	try {
	    ConfigFile config(cfg.prec_conf_file);

	    /*                                          name in file  | default value */
	    qs_str	    = config.read<std::string>	("q",		"0.8"	    );
	    inner_iter_type = config.read<std::string>	("inner_iter",	"fixed"	    );
	    lu_type         = config.read<std::string>	("lu_method",	"exact"	    );
	    elim_order_type = config.read<std::string>	("elim_order",	"original"  );
	    elim_stop_type  = config.read<std::string>	("elim_stop",	"degree"    );
	    coarse_n        = config.read<uint>		("coarse_n",	0	    );
	    nlevels         = config.read<uint>		("max_levels",	15	    );

	    qs		    = new_vector<double>(qs_str);

	    if (inner_iter_type == "fixed") {
		niters_str  = config.read<std::string>	("level_niter", "5"	    );
		niters	    = new_vector<uint>(niters_str);
	    } else if (inner_iter_type == "dynamic") {
		epss_str    = config.read<std::string>	("level_eps",	"1e-2"	    );
		epss	    = new_vector<double>(epss_str);
	    } else
		THROW_EXCEPTION("Unknown inner_iter type: " << inner_iter_type);

	    if (lu_type == "ilut") {
		ilut_p	    = config.read<uint>		("ilut_p",	7	    );
		ilut_tau    = config.read<double>	("ilut_tau",	1e-4	    );
	    } else if (lu_type != "exact")
		THROW_EXCEPTION("Unknown lu type: " << lu_type);

	    if (elim_order_type != "none" &&
		elim_order_type != "original" &&
		elim_order_type != "block" &&
		elim_order_type != "simple_1")
		THROW_EXCEPTION("Unknown elim_order type: " << elim_order_type);

	    degree_max  = 0;
	    level_ratio = 0.0;
	    if (elim_stop_type == "degree")
		degree_max  = config.read<uint>		("degree_max",	3	    );
	    else if (elim_stop_type == "goal") {
		level_ratio = config.read<double>	("level_ratio", 0.25	    );
	    } else
		THROW_EXCEPTION("Unknown elim_stop type: " << elim_stop_type);

	} catch (ConfigFile::file_not_found) {
	    LOG_INFO("Config file \"" << cfg.prec_conf_file << "\" is absent");
	}
    } else { // if (cfg.prec_conf_file)
	inner_iter_type = "fixed";
	lu_type         = "exact";
	elim_order_type = cfg.use_tails ? "original" : "none";
	elim_stop_type  = "degree";
	coarse_n        = cfg.coarse_n;
	nlevels         = cfg.max_levels;

	qs              = cfg.sigmas;
	niters          = cfg.niters;
    }

    if      (elim_order_type == "none")	    order = ORDER_NONE;
    else if (elim_order_type == "original") order = ORDER_ORIGINAL;
    else if (elim_order_type == "block")    order = ORDER_BLOCK;
    else if (elim_order_type == "simple_1") order = ORDER_SIMPLE_1;

    /* Constructing the preconditioner according to parameters */
    levels.resize(nlevels);

    /* Fill in q parameter for each level */
    for (uint l = 0; l < nlevels; l++)
	levels[l].q = (l >= qs.size() ? qs.back() : qs[l]);

    /* Set number of iterations per level */
    levels[0].niter = 1;
    LOG_INFO("The number of iterations on the first level is always 1");
    for (uint l = 1; l < nlevels; l++) {
	if (inner_iter_type == "fixed")
	    levels[l].niter = (l >= niters.size() ? niters.back() : niters[l]);
	else
	    levels[l].eps   = (l >= epss.size()   ? epss.back()   : epss[l]);
    }

    construct_level(0, A);

    levels.resize(nlevels);
    levels[nlevels-1].niter = 0;
    levels[nlevels-1].eps   = 0.0;
    levels[nlevels-1].q	    = 0.0;

#ifdef HAVE_UMFPACK
    Ac_symbolic = Ac_numeric = NULL;
#endif

#ifdef PRINT_NORMS
    /* Initialize stream for dumping norms */
    norm_oss = new std::ostringstream;
    (*norm_oss) << std::scientific;
#endif

    {
	/* Run analysis for A_1 */
	const SkylineMatrix& S = levels[1].A;
	uint nbad = 0;
	float tau = 0.8;
	for (uint i = 0; i < S.size(); i++)
	    if (1 - S.row_sum(i)/S(i,i) > tau)
		nbad++;
	std::cout << "Number of bad nodes for A_1 (tau = " << tau << "): " << nbad << std::endl;
    }
}

MultiSplitPrec::~MultiSplitPrec() {
#ifdef HAVE_UMFPACK
    if (Ac_symbolic && Ac_numeric)
	DirectSolver_free(Ac_symbolic, Ac_numeric);
#endif
#ifdef PRINT_NORMS
    dump_norm_trace();
    delete norm_oss;
#endif
}
