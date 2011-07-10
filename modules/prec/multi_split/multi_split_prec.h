#ifndef __MULTI_SPLIT_PREC_H__
#define __MULTI_SPLIT_PREC_H__

#include "config/config.h"
#include "project/config.h"
#include "modules/prec/prec_base.h"
#include "modules/prec/misc/misc.h"
#include "modules/mesh/mesh.h"

#include <iostream>

#define PRINT_NORMS

class LinkTypeMultiSplit;
class MultiSplitPrec : public PrecBase {
private:
    bool use_tails;

    struct Level {
	uint N, nnz;		    /* Number of nodes and nonzeros elements for the level */
	uint M;                     /* Size of the excluded block (not including diagonal block) */
	uint Md;		    /* Size of the excluded diagonal block */

	double q;		    /* Convergence factor */
	uint niter;		    /* Number of iterations on the level */

	SkylineMatrix A;	    /* Level matrix (for level 0 we use level0_A) */
	CSRMatrix     L;	    /* L factor for the level (N x N) */
	SkylineMatrix U;	    /* U factor for the level (M x N) */

	uvector<uint> map;	    /* Indices map: permuted -> original */
	uvector<uint> rmap;	    /* Indices map: original -> permuted */

	uvector<double> dval;	    /* Reciprocal of the diagonal of the diagonal block */

	mutable
	Vector r, w, x2, u0, F;	    /* Some auxilary vectors for internal iterations */
    };

    uint nlevels;
    std::vector<Level> levels;

    const SkylineMatrix& level0_A;  /* Matrix for level 0 is not kept in levels so that we don't
				       need to make a copy */

    uint coarse_n;
#ifdef HAVE_UMFPACK
    /* UMFPACK factorization of the coarsest level */
    mutable void *Ac_symbolic, *Ac_numeric;
#endif

#ifdef PRINT_NORMS
    mutable std::ostringstream *norm_oss;
    void dump_norm_trace() const;
#endif

    void construct_level(uint level, const SkylineMatrix& A);

    void solve(uint level, const Vector& f, Vector& x) const THROW;

    void solve_U(uint level, const Vector& w, Vector& x) const;
    void solve_L(uint level, const Vector& f, Vector& w) const;
    void solve_diagonal(uint level, const Vector& f, Vector& x) const THROW;

    void construct_permutation(const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
			       uvector<int>& nlinks_in, uvector<int>& nlinks_out,
			       uint& Md, uint& M, uvector<uint>& map) const;

public:
    MultiSplitPrec(const SkylineMatrix& A, const Config& cfg);
    ~MultiSplitPrec();

    void solve(Vector& f, Vector& x) const THROW; /* Wrapper for solve(level,f,x) */
    void graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const;

    friend std::ostream& operator<<(std::ostream& os, const MultiSplitPrec& p);
    friend std::ostream& operator<<(std::ostream& os, const MultiSplitPrec::Level& li);
};

#endif
