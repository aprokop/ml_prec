#include "prec.h"
#include "include/logger.h"

DEFINE_LOGGER("Prec");

Prec::Prec(uint _nlevels, double eps, uint ncheb, double _c, const SkylineMatrix& A) {
    ASSERT(A.size(), "Matrix has size 0");

    nlevels = _nlevels;
    levels.resize(nlevels);
    c = _c;

    // alpha and beta for each level
    for (uint level = 0; level < nlevels-1; level++) {
	levels[level].alpha = 1.;
	levels[level].beta  = eps;
	levels[level].ncheb = ncheb;
    }

    construct_level(0, A);
}

void Prec::construct_level(uint level, const SkylineMatrix& A) {
    Level& li = levels[level];
    uint N = li.N = A.size();

    // ASSERT(A.is_symmetric(), "Level: " << level << ", A is not symmetric");

    if (level < nlevels-1) {
	SkylineMatrix& nA = levels[level+1].A;
	std::vector<uint>& tr = li.tr;
	std::vector<uint>& dtr = li.dtr;

	// reverse translation
	std::vector<uint> revtr(N, -1);

	nA.ia.push_back(0);
	for (uint i = 0; i < N; i++) {
	    // dind corresponds to the position of diagonal
	    uint dind = nA.a.size(); 

	    nA.ja.push_back(i);
	    nA.a.push_back(c);

	    double v;
	    for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++) 
		if (1 + 12*(-A.a[j]) / c > li.beta) {
		    // if the link stays, scale it
		    nA.ja.push_back(A.ja[j]);
		    v = A.a[j] / li.beta;
		    nA.a.push_back(v);
		    nA.a[dind] += -v;
		}
	    
	    uint dia = nA.a.size() - dind;
	    if (dia > 1) {
		nA.ia.push_back(nA.a.size());

		tr.push_back(i);
		revtr[i] = nA.ia.size()-2;
	    } else {
		// all links for this point are gone
		nA.ja.pop_back();
		nA.a.pop_back();

		// add point to the diagonal vector
		dtr.push_back(i);
	    }
	}

	// change nA.ja to use local indices
	for (uint j = 0; j < nA.ja.size(); j++) {
	    ASSERT(revtr[nA.ja[j]] != uint(-1), "Trying to invert wrong index: j = " << j << ", nA.ja[j] = " << nA.ja[j]);
	    nA.ja[j] = revtr[nA.ja[j]];
	}
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
	}
    }

    if (level == nlevels-1) {
	// this parameters do not matter, it is not used; just want to set them to smth
	Level& lc = levels[nlevels-1];
	lc.ncheb = 0; 
	lc.alpha = lc.beta = 0.;

	// these DO matter
	lc.lmin = lc.lmax = 1;

	levels[nlevels-2].ncheb = 0;

	// chebyshev info
	for (int l = nlevels-2; l >= 0; l--) {
	    Level& li = levels[l];
	    Level& ln = levels[l+1];

	    if (li.ncheb) {
		double cs = cheb((ln.lmax + ln.lmin)/(ln.lmax - ln.lmin), li.ncheb);
		li.lmin = li.alpha*(1 - 1/cs);
		li.lmax = li.beta *(1 + 1/cs);
	    } else {
		li.lmin = li.alpha*ln.lmin;
		li.lmax = li.beta *ln.lmax;
	    }
	}
    }
}
