#include "prec.h"
#include "include/logger.h"

#include <map>

DEFINE_LOGGER("Prec");

Prec::Prec(double eps, uint _ncheb, const SkylineMatrix& A, const Mesh& _mesh) : mesh(_mesh) {
    ASSERT(A.size(), "Matrix has size 0");

    galpha = 1.;
    gbeta  = eps;
    ncheb = _ncheb;

    // reserve
    nlevels = 20;
    levels.resize(nlevels);
    levels[0].A = A;
    construct_level(0, A);
}

class LinkType {
private:
    std::vector<std::map<uint,char> > vec;

    void check(uint i, uint j) const {
	ASSERT(i != j, "i == j == " << i);
	ASSERT(i < vec.size() && j < vec.size(), "i = " << i << ", j = " << j << ", vec.size() == " << vec.size());
    }

public:
    LinkType(uint n) : vec(n) { }

    char operator()(uint i, uint j) const {
	check(i,j);

	uint i0, i1;
	if (i < j) { i0 = i; i1 = j; }
	else       { i0 = j; i1 = i; }

	std::map<uint,char>::const_iterator it = vec[i0].find(i1);
	if (it != vec[i0].end())
	    return it->second;
	return 0;
    }

    char& operator()(uint i, uint j) {
	check(i,j);
	if (i < j)
	    return vec[i][j];
	return vec[j][i];
    }
};

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
    std::multimap<double,uint> rmap;
    for (uint i = 0; i < N; i++) {
	rmap.clear();

	// aux[i] == c[i]
	aux[i] = A.a[A.ia[i]];
	for (uint j = A.ia[i]+1; j < A.ia[i+1]; j++) {
	    rmap.insert(std::pair<double,uint>(-A.a[j], A.ja[j]));
	    aux[i] += A.a[j];
	}

	nlinks[i] = A.ia[i+1] - A.ia[i] - 1;

	double s = 0.;
	for (std::multimap<double,uint>::const_iterator it = rmap.begin(); it != rmap.end(); it++) {
	    s += 2*it->first / (aux[i]*(gbeta-1));
	    if (s <= 1) {
		uint i0, i1;
		ltype(i,it->second)++;

		if (i > it->second && ltype(i,it->second) == 2) {
		    nlinks[i]--;
		    nlinks[it->second]--;
		}
	    } else {
		break;
	    }
	}
    }

    // ===============  STEP 2 : link removing : not using c  ===============
    // it changes spectrum on the finest level to at least [1,4], but we don't use that estimate
#if 1
    if (level == 0) {
	uint removed = 0;

	uint iis[4];
	uint nx = mesh.nx, ny = mesh.ny, nz = mesh.nz;
#if 0
	for (uint k = 0; k < nz; k++)
	    for (uint j = 0; j < ny-1; j++)
		for (uint i = j&1; i < nx-1; i += 2) {
		    // xy plane
		    iis[0] = mesh.index(  i,  j,k);
		    iis[1] = mesh.index(i+1,  j,k);
		    iis[2] = mesh.index(i+1,j+1,k);
		    iis[3] = mesh.index(  i,j+1,k);
#else
	for (uint i = 0; i < nx; i++) 
	    for (uint k = 0; k < nz-1; k++)
		for (uint j = k&1; j < ny-1; j += 2) {
		    // yz plane
		    iis[0] = mesh.index(i,  j,  k);
		    iis[1] = mesh.index(i,j+1,  k);
		    iis[2] = mesh.index(i,j+1,k+1);
		    iis[3] = mesh.index(i,  j,k+1);
#endif

		    if (ltype(iis[0],iis[1]) == 2 ||
			ltype(iis[1],iis[2]) == 2 ||
			ltype(iis[2],iis[3]) == 2 ||
			ltype(iis[3],iis[0]) == 2)
			continue;

		    uint i0, i1;
#if 0
		    double min = -A.get(iis[0],iis[1]), t;
		    if ((t = -A.get(iis[1],iis[2])) < min) { min = t; i0 = 1; i1 = 2; }
		    if ((t = -A.get(iis[2],iis[3])) < min) { min = t; i0 = 2; i1 = 3; }
		    if ((t = -A.get(iis[3],iis[0])) < min) { min = t; i0 = 3; i1 = 0; }
#else
		    rmap.clear();
		    // LOG_DEBUG("iis: " << iis[0] << ", " << iis[1] << ", " << iis[2] << ", " << iis[3]);
		    // LOG_DEBUG("As : " << -A.get(iis[0],iis[1]) << ", " << -A.get(iis[1],iis[2]) << ", " << 
			      // -A.get(iis[2],iis[3]) << ", " << -A.get(iis[3],iis[0]));

		    rmap.insert(std::pair<double,uint>(-A.get(iis[0],iis[1]), 0));
		    rmap.insert(std::pair<double,uint>(-A.get(iis[1],iis[2]), 1));
		    rmap.insert(std::pair<double,uint>(-A.get(iis[2],iis[3]), 2));
		    rmap.insert(std::pair<double,uint>(-A.get(iis[3],iis[0]), 3));

		    std::multimap<double,uint>::const_iterator it = rmap.begin();
		    double min = it->first, smin = (++it)->first;
		    // LOG_DEBUG("(1) " << min << ", (2) " << smin);
		    if (smin < 3./(gbeta-1)*min)
			continue;
		    i0 = rmap.begin()->second;
#endif
		    i1 = (i0+1)&3;

		    ltype(iis[i0],iis[i1]) = 2;
		    nlinks[iis[i0]]--;
		    nlinks[iis[i1]]--;

		    removed++;
		}
	std::cout << "removed: " << removed << std::endl;
	ASSERT(false, "want to abort");
    }
#endif
    
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
