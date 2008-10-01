#include "rel_prec.h"
#include "include/logger.h"
#include "include/tools.h"

#include <numeric>
#include <map>

DEFINE_LOGGER("Prec");

RelPrec::RelPrec(const SkylineMatrix& A, uint _niter, double _gamma, 
		 const std::vector<double>& _sigmas, const Mesh& _mesh) : mesh(_mesh) {
    ASSERT(A.size(), "Matrix has size 0");
    ASSERT(_gamma >= 1, "Wrong gamma: " << _gamma);
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

void RelPrec::construct_level(uint level, const SkylineMatrix& A) {
    ASSERT(level < nlevels, "Incorrect level: " << level << " (" << nlevels << ")");
    Level& li = levels[level];

    // ASSERT(A.is_symmetric(), "Level: " << level << ", A is not symmetric");

    li.N   = A.size();
    li.nnz = A.ja.size();

    uint N = li.N;
    std::vector<double>& aux = li.aux;
    aux.resize(N);

    std::vector<int> nlinks(N);
    LinkType ltype(N);

    // ===============  STEP 1 : link removing : using c  ===============
    // calculate c
    for (uint i = 0; i < N; i++) 
	aux[i] = std::accumulate(A.a.begin() + A.ia[i], A.a.begin() + A.ia[i+1], 0.);

    typedef std::multimap<double,uint> map_type;
    map_type rmap;
    for (uint i = 0; i < N; i++) {
	rmap.clear();
	for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++) 
	    rmap.insert(std::pair<double,uint>(-A.a[j], A.ja[j]));
	nlinks[i] = A.ia[i+1] - A.ia[i] - 1;

	double c = aux[i];

	double s = 0., beta;
	map_type::const_iterator it;
	for (it = rmap.begin(); it != rmap.end(); it++) {
	    // s + 2a/(c(sigma-1)) <= 1
	    std::vector<double>::const_iterator lb = std::lower_bound(sigmas.begin(), sigmas.end(), 
								      1 + 2*it->first / (c*(1-s)));
	    if (lb != sigmas.end()) {
		double sigma = *lb;
		beta = 2*it->first / (c*(sigma-1));
		s += beta;
		uint i0, i1;

		if (i < it->second) {
		    // scale this part by gamma in case 
		    // link won't be removed later
		    aux[i] += (gamma-1)*beta*c;
		    ltype(i,it->second)++;
		}
		if (i > it->second) {
		    if (ltype(i,it->second) == 1) {
			// remove link
			ltype(i,it->second) = 2;

			nlinks[i]--;
			nlinks[it->second]--;
			aux[i]          += (sigma - 1)*beta*c;
			aux[it->second] += (sigma - gamma)*beta*c;
		    } else {
			// we don't delete link
			aux[i] += (gamma-1)*beta*c;
		    }
		}
	    } else {
		aux[i] += (gamma - 1)*(1-s)*c;
		break;
	    }
	}
    }
    
    // ===============  STEP 2 : tail removing ===============
    std::vector<Tail>& tails = li.tails;
    tails.clear(); // just in case

    uint i0 = -1, i1 = -1;
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

    std::vector<uint>& tr  = li.tr;
    std::vector<uint>& dtr = li.dtr;

    // reverse translation
    std::vector<uint> revtr(N, -1);

    // ===============  STEP 3 : construct matrix  ===============
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

	li.x0.resize(n);
	li.x1.resize(n);
	li.f1.resize(n);

	construct_level(level+1, nA);
    } else {
	LOG_INFO("Decreasing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level+1;
	levels.resize(nlevels);
    }
}
