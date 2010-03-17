#ifndef __MULTI_SPLIT_PREC_H__
#define __MULTI_SPLIT_PREC_H__

#include "config/config.h"
#include "modules/common/common.h"
#include "modules/prec/prec_base.h"
#include "modules/mesh/mesh.h"

#include <iostream>

class MultiSplitPrec : public PrecBase {
private:
    bool use_tails;

    struct Level {
	uint N, nnz;		    /* Number of nodes and nonzeros elements for the level */
	double q;		    /* Convergence factor */

	uint niter;		    /* Number of iterations on the level */

	SkylineMatrix A;	    /* Level matrix (for level 0 we use level0_A) */

	mutable
	Vector x1, f1, u0, r;	    /* Some auxilary vectors for internal iterations */

	uvector<uint> tr;	    /* C-indices */
	uvector<uint> dtr;	    /* F-indices */

	uvector<double> aux;	    /* Auxilary array, corresponds to value c of the node; also
				       used in tail elimination */
    };
    uint nlevels;
    std::vector<Level> levels;

    const SkylineMatrix& level0_A;  /* Matrix for level 0 is not kept in levels so that we don't
				       need to make a copy */

    void construct_level(uint level, const SkylineMatrix& A);

    void solve(uint level, Vector& f, Vector& x) const THROW;

public:
    MultiSplitPrec(const SkylineMatrix& A, const Config& cfg);
    ~MultiSplitPrec() { }

    void solve(Vector& f, Vector& x) THROW; /* Wrapper for solve(level,f,x) */
    void graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const;

    friend std::ostream& operator<<(std::ostream& os, const MultiSplitPrec& p);
    friend std::ostream& operator<<(std::ostream& os, const MultiSplitPrec::Level& li);
};

#endif
