#include "prec.h"
#include "include/logger.h"

#include <set>

DEFINE_LOGGER("Prec");

Prec::Prec(double eps, uint _ncheb, double _c, const SkylineMatrix& A) {
    ASSERT(A.size(), "Matrix has size 0");

    galpha = 1.;
    gbeta  = eps;
    c = _c;
    ncheb = _ncheb;

    // reserve
    nlevels = 20;
    levels.resize(nlevels);
    construct_level(0, A);
}

void Prec::construct_level(uint level, const SkylineMatrix& A) {
    ASSERT(level < nlevels, "Incorrect level: " << level << " (" << nlevels << ")");
    Level& li = levels[level];

    // ASSERT(A.is_symmetric(), "Level: " << level << ", A is not symmetric");

    li.N   = A.size();
    li.nnz = A.ja.size();
    uint N = li.N;

    std::vector<uint> nlinks(N);
    std::vector<double> _c(N);

    std::vector<std::map<uint,char> > vec(N);
    std::multimap<double,uint> rmap;
    for (uint i = 0; i < N; i++) {
	rmap.clear();

	_c[i] = A.a[A.ia[i]];
	for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++) {
	    rmap.insert(std::pair<double,uint>(-A.a[j], A.ja[j]));
	    _c[i] += A.a[j];
	}

	nlinks[i] = A.ia[i+1] - A.ia[i] - 1;

	double s = 0.;
	for (std::multimap<double,uint>::const_iterator it = rmap.begin(); it != rmap.end(); it++) {
	    s += 2*it->first / (_c[i]*(gbeta-1));
	    if (s <= 1) {
		uint i0, i1;
		if (i < it->second) { i0 = i; i1 = it->second; }
		else		    { i0 = it->second; i1 = i; }
		vec[i0][i1]++;

		if (i > it->second && vec[i0][i1] == 2) {
		    nlinks[i]--;
		    nlinks[it->second]--;
		}
	    } else {
		break;
	    }
	}
    }
    
    for (uint i = 0; i < N; i++) 
	if (nlinks[i] == 1) {
	    uint i0 = i, i1;
	    do {
		for (uint j = A.ia[i0]+1; j < A.ia[i0+1]; j++) {
		    uint jj = A.ja[j];
		    if (i0 < jj && vec[i0][jj] != 2 || i0 > jj && vec[jj][i0] != 2) {
			i1 = jj;
			break;
		    }
		}
		nlinks[i0]--;
		nlinks[i1]--;
		if (i0 < i1) vec[i0][i1] = 2;
		else         vec[i1][i0] = 2;
		i0 = i1;
	    } while (nlinks[i0] == 1);
	}

    std::vector<uint>& tr  = li.tr;
    std::vector<uint>& dtr = li.dtr;

    // reverse translation
    std::vector<uint> revtr(N, -1);

    SkylineMatrix& nA = levels[level+1].A;
    nA.ia.push_back(0);
    for (uint i = 0; i < N; i++) 
	if (nlinks[i]) {
	    // this node goes to the next level

	    // dind corresponds to the position of diagonal
	    uint dind = nA.a.size(); 

	    nA.ja.push_back(i);
	    nA.a.push_back(c);

	    for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++) {
		uint jj = A.ja[j];
		// dynamic c
		if (i < jj && vec[i][jj] != 2 || i > jj && vec[jj][i] != 2) {
		    // link stays, scale it
		    nA.ja.push_back(A.ja[j]);
		    double v = A.a[j] / gbeta;
		    nA.a.push_back(v);
		    nA.a[dind] += -v;
		}
	    }

	    nA.ia.push_back(nA.a.size());

	    tr.push_back(i);
	    revtr[i] = nA.ia.size()-2;
	} else {
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

    // add local indices numbers to tails

    uint n = tr.size();
    if (n) {
	nA.nrow = nA.ncol = n;
	revtr.clear();

	li.x1.resize(n);
	li.u0.resize(n);
	li.u1.resize(n);
	li.f1.resize(n);

	construct_level(level+1, nA);
    } else {
	LOG_INFO("Decreasing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level+1;
	levels.resize(nlevels);

	Level& lc = levels[nlevels-1];
	lc.lmin = 1;
	lc.lmax = gbeta;

	// chebyshev info
	for (int l = nlevels-2; l >= 0; l--) {
	    Level& li = levels[l];
	    Level& ln = levels[l+1];

	    double cs = cheb((ln.lmax + ln.lmin)/(ln.lmax - ln.lmin), ncheb);
	    li.lmin = galpha*(1 - 1/cs);
	    li.lmax = gbeta *(1 + 1/cs);
	}
    }
}
