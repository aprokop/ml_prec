#include "prec.h"
#include "include/logger.h"

DEFINE_LOGGER("Prec");

Prec::Prec(uint nlevels, double eps, uint ncheb, const SkylineMatrix& A) {
    ASSERT(A.size(), "Matrix has size 0");
    levels.resize(nlevels);

    // alpha and beta for each level
    for (int level = 0; level < nlevels-1; level++) {
	levels[level].alpha = 1.;
	levels[level].beta  = eps;
	levels[level].ncheb = ncheb;
    }
    // this parameters do not matter, it is not used; just want to set them to smth
    Level& lc = levels[nlevels-1];
    lc.ncheb = 0; 
    lc.alpha = lc.beta = 0.;

    construct_level(0, A);
}

void Prec::construct_level(uint i, const SkylineMatrix& A) {
    Level& li = levels[i];
    uint N = li.N = A.size();

    if (i < nlevel-1) {
	SkylineMatrix& nA = levels[i+1].A;

	// reverse translation
	std::vector<int> revtr(N, -1);

	for (uint i = 0; i < N; i++) {
	    nA.ja.push_back(i);
	    nA.a.push_back(c);

	    // ind corresponds to the position of diagonal
	    uint dind = nA.a.size()-1; 

	    for (int j = A.ia[i]+1; j < A.ia[i+1]; j++) {
		if (1 + 12*(-A.a[j]) / c > eps) {
		    // if the link stays, scale it
		    nA.ja.push_back(A.ja[j]);
		    nA.a.push_back(A.a[j] / li.beta);
		    nA.a[dind] += (-A.a[j]) / li.beta;
		}
	    }
	    uint dia = nA.a.size() - dind;
	    if (dia != 1) {
		uint ni = nA.ia.size() - 1;
		nA.ia.push_back(nA.ia[ni] + dia);
		tr.push_back(i);

		revtr[i] = ni;
	    } else {
		// all links for this point are gone
		nA.ja.pop_back();
		nA.a.pop_back();

		// add point to the diagonal vector
		dtr.push_back(i);
	    }
	}
	// change nA.ja to use local indices
	for (int j = 0; j < nA.ja.size(); j++) {
	    nA.ja[j] = revtr[nA.ja[j]];
	}
	nA.nrow = nA.ncol = tr.size();
	revtr.clear();

	x1.resize(tr.size());
	f1.resize(tr.size());

	construct_level(i+1, nA);
    } else {
	// chebyshev info
	for (int level = nlevels-2; level >= 0; level--) {
	    Level& li = levels[level];
	    Level& ln = levels[level+1];

	    double cs = cheb((ln.lmax + ln.lmin)/(ln.lmax - ln.lmin), li.ncheb);
	    li.lmin = li.alpha*(1 - 1/cs);
	    li.lmax = li.beta *(1 + 1/cs);
	}
    }
}

void Prec::solve(const Vector& f, Vector& x) const THROW {
    solve(f, x, 0);
}

void Prec::solve(const Vector& f, Vector& x, uint level) const THROW {
    int N = levels[level].N;
    ASSERT(f.size() == N && x.size() == N, "Wrong dimension: N = " << n << ", f = " << f.size() << ", x = " << x.size());

    if (level < nlevels-1) {
	uint n = levels[level+1].N;

	for (int i = 0; i < n; i++) 
	    f1[i] = f[tr[i]];

	solve(x1, f1, level+1);

	for (int i = 0; i < n; i++)
	    x[tr[i]] = x1[i];

    } else {
	// for last level assert for now that we have only diagonal matrix
	ASSERT(dtr.size() == N, "Now we must have a diagonal on the coarsest level");
    }

    // solve diagonal part
    for (int i = 0; i < dtr.size(); i++)
	x[dtr[i]] = f[dtr[i]] / c;
}

std::ostream& operator<<(std::ostream& os, const Prec& p) {
    os << "nlevels = " << p.nlevels;
    for (int level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	const Prec::Level& li = p.levels[level];
	os << "N = " << li.N << std::endl;
	if (level != p.nlevels-1)
	    os << "Ncheb = " << li.ncheb << ", ";
	os << "[lmin, lmax] = [" << li.lmin << "," << li.lmax << "]" << std::endl;
	os << "alpha = " << li.alpha << ", beta = " << li.beta << std::endl;
#if 1
	if (level < nlevels-1) {
	    os << "tr: " << li.tr;
	    os << "dtr: " << li.dtr;
	}
#endif
#if 0
	if (level) 
	    os << "A: " << li.A;
#endif
    }
    if (p.levels[p.nlevels-1].ncheb)
	std::cout << "Coarse ncheb = " << p.levels[p.nlevels-1].ncheb << std::endl;
    return os;
} 
