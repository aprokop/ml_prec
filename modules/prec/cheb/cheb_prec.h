#ifndef __PREC_H__
#define __PREC_H__

#include "config/config.h"
#include "cheb_misc.h"
#include "modules/prec/misc/misc.h"
#include "modules/common/common.h"
#include "modules/prec/prec_base.h"
#include "modules/mesh/mesh.h"
#include <iostream>

class Prec : public PrecBase {
private:
    bool use_tails;

    struct Level {
	uint N, nnz;		    /* Number of nodes and nonzeros elements for the level */

	double alpha, beta;	    /* Spectral constants used for the level */
	double lmin, lmax;	    /* Spectral boundaries for the level */
	uint ncheb;		    /* Number of Chebyshev iterations on the level */

	SkylineMatrix A;	    /* Level matrix (for level 0 we use level0_A) */

	mutable
	Vector x1, f1, u0, u1, tmp; /* Some auxilary vectors for Chebyshev iterations */

	uvector<uint> tr;	    /* C-indices */ 
	uvector<uint> dtr;	    /* F-indices */
	std::vector<Tail> tails;    /* Tridiagonal matrices */

	uvector<double> aux;	    /* Auxilary array, corresponds to value c of the node; also 
				       used in tail elimination */
    };
    uint nlevels;
    std::vector<Level> levels;

    const SkylineMatrix& level0_A;  /* Matrix for level 0 is not kept in levels so that we don't
				       need to make a copy */

    /* 
     * This function tries to increase s by amount necessary to mark v as removable.
     * If it can not be done (i.e. there is not enough value of c for that then s
     * is not changed
     * NOTE: if matrix contains positive offdiagonal elements (as in the
     * case when we take Schur complement to u and p in mixed hybrid FEM)
     * removing negative offdiagonal elements may result in matrix which is
     * not positive definite. Sometimes we can assume that if c is big 
     * enough new matrix would still be positive definite
     */
    bool to_remove(double c, double v, double beta, double& s) {
	double t = -2*v / (c*(beta-1));
#if 0
	/* There can be several possibilities */
	if (v < 0 && beta > 1) {
	    /* We are allowed to touch negative off-diagonal elements */
	    t = -2*v / (c*(beta-1));
	} else if (v > 0 && alpha < 1) {
	    /* We are allowed to touch positive off-diagonal elements */
	    t =  2*v / (c*(1-alpha));
	} else {
	    // LOG_WARN("Zero offdiagonal element");
	} 
#endif

	if (s+t > 1)
	    return false;

	s += t;
	return true;
    }

    void optimize_storage();
    void construct_level(uint level, const SkylineMatrix& A);

    void solve(uint level, Vector& f, Vector& x) const THROW;

public:
    Prec(const SkylineMatrix& A, const Config& cfg);
    ~Prec() { }

    void solve(Vector& f, Vector& x) THROW; /* Wrapper for solve(level,f,x) */
    void graph_planes(const std::string& filename, uint level, char plane) const;

    double lmin() const { return levels[0].lmin; }
    double lmax() const { return levels[0].lmax; }

    friend std::ostream& operator<<(std::ostream& os, const Prec& p);
    friend std::ostream& operator<<(std::ostream& os, const Prec::Level& li);
};

#endif
