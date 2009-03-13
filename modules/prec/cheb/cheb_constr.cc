#include "include/logger.h"
#include "cheb_prec.h"
#include "modules/common/common.h"

#include <numeric>
#include <map>

DEFINE_LOGGER("Prec");

Prec::Prec(double eps, uint _ncheb, const SkylineMatrix& A, const MeshBase& _mesh) : mesh(_mesh) {
    ASSERT(A.size(), "Matrix has size 0");

    switch (2) {
	case 1:
	    galpha = 1;
	    gbeta  = eps;
	    break;
	case 2:
	    gbeta = .8*eps;
	    galpha = gbeta/eps;
	    break;
    }
    ncheb = _ncheb;

    // reserve
    nlevels = 20;
    levels.resize(nlevels);
    levels[0].A = A;
    construct_level(0, A);
}

// insertion sort of sorted wrt to abs values of a
static void psort(const double * a, uint n, std::vector<uint>& sorted) {
    uint i;
    int j;
    for (uint i = 0; i < n; i++)
	sorted[i] = i;

    uint c;
    double v;
    for (i = 1; i < n; i++) {
	v = a[sorted[i]];
	c = sorted[i];
	for (j = i-1; j >= 0 && a[sorted[j]] > v; j--)
	    sorted[j+1] = sorted[j];
	sorted[j+1] = c;
    }
}

void Prec::construct_level(uint level, const SkylineMatrix& A) {
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
    int rstart, rend;
    std::vector<uint> sorted;
    for (uint i = 0; i < N; i++) {
	rstart = A.ia[i];   // row start
	rend   = A.ia[i+1]; // row end

	// aux[i] == c[i]
	aux[i] = 0;
	for (int j = rstart; j < rend; j++)
	    aux[i] += A.a[j];
	nlinks[i] = rend - rstart - 1;
	sorted.resize(nlinks[i]);

	// sort off-diagonal elements wrt to their abs values
	psort(&A.a[rstart+1], nlinks[i], sorted);

	double s = 0.;
	for (int j = 0; j < nlinks[i]; j++) {
	    double a = A.a[rstart+1 + sorted[j]];
	    if (a < 0) {
		// for now we always try to eliminate negative components
		s += -2*a / (aux[i]*(gbeta-1));
	    } else if (galpha != 1) {
		// if we are allowed to change offdiagonal positive
		s +=  2*a / (aux[i]*(1-galpha));
	    } else {
		// offdiagonal element is positive but we are not allowed to change it
		continue;
	    }

	    if (s <= 1) {
		uint i0 = i, i1 = A.ja[rstart+1 + sorted[j]];
		ltype(i0, i1)++;

		if (i0 > i1 && ltype(i0,i1) == 2) {
		    nlinks[i0]--;
		    nlinks[i1]--;
		}
	    } else {
		break;
	    }
	}
    }

    // ===============  STEP 3 : tail removing ===============
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
			v = A.a[j] / gbeta;
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
		// dynamic c
		if (ltype(i,jj) != 2) {
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

	li.x1.resize(n);
	li.f1.resize(n);
	li.u0.resize(n);
	li.u1.resize(n);

	construct_level(level+1, nA);
    } else {
	LOG_INFO("Decreasing number of levels: " << nlevels << " -> " << level+1);
	nlevels = level+1;
	levels.resize(nlevels);

	Level& lc = levels[nlevels-1];
	lc.lmin = galpha;
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
