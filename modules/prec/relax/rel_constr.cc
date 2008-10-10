#include "rel_prec.h"
#include "include/logger.h"
#include "include/tools.h"
#include "include/time.h"

#include <cstring>
#include <algorithm>
#include <numeric>
#include <map>
#include <vector>

DEFINE_LOGGER("RelPrec");

RelPrec::RelPrec(const SkylineMatrix& A, uint _niter, double _gamma, 
		 const std::vector<double>& _sigmas, const MeshBase& _mesh) : mesh(_mesh) {
    ASSERT(A.size(), "Matrix has size 0");
    ASSERT(_gamma > 1, "Wrong gamma: " << _gamma);
    ASSERT(_sigmas.size(), "Sigmas is empty");
    for (uint i = 0; i < _sigmas.size(); i++)
	ASSERT(_sigmas[i] > 1, "Wrong sigmas[" << i << "]: " << _sigmas[i]);

    niter  = _niter;
    gamma  = _gamma;
    sigmas = _sigmas;

    // reserve
    nlevels = 20;
    levels.resize(nlevels);
    levels[0].A = A;
    construct_level(0, A);
}

struct Link {
    double val;
    uint i0, i1;

    Link(double v, uint i, uint j) : val(v), i0(i), i1(j) { }
    friend bool operator<(const Link& l0, const Link& l1) {
	return l0.val < l1.val;
    }
};

void RelPrec::construct_level(uint level, const SkylineMatrix& A) {
    ASSERT(level < nlevels, "Incorrect level: " << level << " (" << nlevels << ")");
    Level& li = levels[level];

    // ASSERT(A.is_symmetric(), "Level: " << level << ", A is not symmetric");
    TIME_INIT();

    li.N   = A.size();
    li.nnz = A.ja.size();

    uint N = li.N;
    Vector& aux = li.aux;
    aux.resize(N);
    Vector remc(N);

    // ! nlinks must be <int>
    std::vector<int> nlinks(N);
    std::vector<Link> links;
    links.reserve(N);
    LinkType ltype(N);

    // ===============  STEP 1 : construct links  ===============
    TIME_START();
    std::vector<uint>::const_iterator start, end, it;
    for (uint i = 0; i < N; i++) {
	nlinks[i] = A.ia[i+1] - A.ia[i] - 1;
	for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
	    remc[i] += A.a[j];

	start = A.ja.begin() + A.ia[i]+1;
	end   = A.ja.begin() + A.ia[i+1];
	it    = std::lower_bound(start, end, i);

	uint jstart = std::distance(A.ja.begin(), it);
	for (uint j = jstart; j < A.ia[i+1]; j++) 
	    links.push_back(Link(-A.a[j], i, A.ja[j]));
    }
    std::sort(links.begin(), links.end());
    LOG_DEBUG(TIME_INFO("# " << level << ": constructing links array"));

    // ===============  STEP 2 : remove links : using c  ===============
    uint i0 = -1, i1 = -1;
    TIME_START();
    for (std::vector<Link>::const_iterator it = links.begin(); it != links.end(); it++) {
	i0 = it->i0;
	i1 = it->i1;
	double c1 = remc[i0]; 
	double c2 = remc[i1];
	if (is_equal(c1, 0.) || is_equal(c2, 0.))
	    continue;
	// calculate what would be the minimum sigma to remove the link
	std::vector<double>::const_iterator lb = std::lower_bound(sigmas.begin(), sigmas.end(), 
								  1 + (c1 + c2)/(c1*c2)*it->val);

	if (lb == sigmas.end()) 
	    continue;

	// remove link with sigma = *lb
	double sigma = *lb;

	// remove link meta data
	nlinks[i0]--;
	nlinks[i1]--;
	ltype(i0, i1) = 2;

	double pc = 2*it->val/(sigma-1), pc1, pc2;
	// note: here we can play some games with
	// distributions of pc1 and pc2 between 
	// c1 and c2
	if (c1 < pc) {
	    pc1 = c1;
	    pc2 = 1/(2/pc - 1/pc1);
	} else if (c2 < pc) {
	    pc2 = c2;
	    pc1 = 1/(2/pc - 1/pc2);
	} else {
	    pc1 = pc2 = pc;
	}

	// subtract c part corresponding to link
	remc[i0] -= pc1;
	remc[i1] -= pc2;

	// update aux array
	aux[i0] += sigma*pc1;
	aux[i1] += sigma*pc2;
    }
    LOG_DEBUG(TIME_INFO("# " << level << ": removing links"));

    // multiply the rest of c by gamma
    daxpy(gamma, remc, aux); 

    // ===============  STEP 2 : tail removing ===============
    TIME_START();
    std::vector<Tail>& tails = li.tails;
    tails.clear(); // just in case

    double v = 0;
#if 1
    for (uint i = 0; i < N; i++) 
	if (nlinks[i] == 1) {
	    Tail tail;

	    i0 = i;
	    do {
		TailNode tn;
		tn.index = i0;
		tn.a3 = -v; // old v

		// find the remaining link
		uint j;
		for (j = A.ia[i0]+1; j < A.ia[i0+1]; j++) {
		    uint jj = A.ja[j];
		    if (ltype(i0,jj) != 2) {
			i1 = jj;
			v = A.a[j];
			break;
		    }
		}
		ASSERT(j < A.ia[i0+1], "??");

		tn.a2  = 1/(aux[i0] - v);
		tn.a1  = -v*tn.a2; // new v
		tn.a3 *= tn.a2;

		nlinks[i0] = -1;
		nlinks[i1]--;

		aux[i1] += aux[i0]*(-v) / (aux[i0] - v);

		ltype(i0,i1) = 2;
		i0 = i1;

		tail.push_back(tn);
	    } while (nlinks[i0] == 1);

	    TailNode tn;
	    tn.index = i0;

	    if (nlinks[i0] == 0) {
		// end node in fully tridiagonal matrix
		tn.a2 = 1/aux[i0];
		nlinks[i0] = -1;
		tail.end_type = 'f';
	    } else {
		tn.a2 = 1.;
		tail.end_type = 'l';
	    }
	    tn.a3 = -v*tn.a2;

	    tail.push_back(tn);

	    tails.push_back(tail);
	}
#endif
    LOG_DEBUG(TIME_INFO("# " << level << ": removing tails"));

    std::vector<uint>& tr  = li.tr;
    std::vector<uint>& dtr = li.dtr;

    // reverse translation
    std::vector<uint> revtr(N, -1);

    // ===============  STEP 3 : construct matrix  ===============
    TIME_START();
    SkylineMatrix& nA = levels[level+1].A;
    nA.ia.push_back(0);
    for (uint i = 0; i < N; i++) 
	if (nlinks[i] > 0) {
	    // the node goes to the next level

	    // dind corresponds to the position of diagonal
	    uint dind = nA.a.size(); 

	    nA.ja.push_back(i);
	    nA.a.push_back(aux[i]);

	    for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++) {
		uint jj = A.ja[j];
		if (ltype(i,jj) != 2) {
		    // link stays
		    nA.ja.push_back(A.ja[j]);
		    double v = A.a[j];
		    nA.a.push_back(v);
		    nA.a[dind] += -v;
		}
	    }

	    nA.ia.push_back(nA.a.size());

	    tr.push_back(i);
	    revtr[i] = nA.ia.size()-2;
	} else if (nlinks[i] != -1) {
	    // all links for this point are gone
	    // add point to the diagonal vector
	    dtr.push_back(i);
	}

    // change nA.ja to use local indices
    for (uint j = 0; j < nA.ja.size(); j++) {
	ASSERT(revtr[nA.ja[j]] != uint(-1), 
	       "Trying to invert wrong index: j = " << j << ", nA.ja[j] = " << nA.ja[j]);
	nA.ja[j] = revtr[nA.ja[j]];
    }
    LOG_DEBUG(TIME_INFO("# " << level << ": constructing matrix"));

    // check whether the ends are really local
    for (uint i = 0; i < tails.size(); i++) {
	Tail& tail = tails[i];
	if (tail.end_type == 'l' && revtr[tail.back().index] == uint(-1))
	    tail.end_type = 't';
    }

    uint n = tr.size();
    if (n) {
	nA.nrow = nA.ncol = n;
	revtr.clear();

	li.tmp0.resize(n);
	li.tmp1.resize(n);
	li.x1.resize(n);
	li.f1.resize(n);

	construct_level(level+1, nA);
    } else {
	LOG_INFO("Decreasing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level+1;
	levels.resize(nlevels);
    }
}
