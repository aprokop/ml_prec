#ifndef __MULTI_SPLIT_PREC_H__
#define __MULTI_SPLIT_PREC_H__

#include "config/config.h"
#include "project/config.h"
#include "modules/prec/prec_base.h"
#include "modules/prec/misc/misc.h"
#include "modules/mesh/mesh.h"

#include <iostream>

// #define PRINT_NORMS

class LinkTypeMultiSplit;
class MultiSplitPrec : public PrecBase {
private:
    bool use_tails;

    struct Level {
        // Config variables
        double	q;		        // convergence factor
        uint	niter;		    // number of iterations on the level
        double	eps;		    // epsilon for residual norm stop criteria

        // Info variables
        uint N, nnz;		    // number of nodes and nonzeros elements for the level
        uint M;                 // size of the excluded block (not including diagonal block)
        uint Md;		        // size of the excluded diagonal block

        // Data variables
        SkylineMatrix	A;	    // level matrix (for level 0 we use level0_A)
        CSRMatrix	    L;	    // L factor for the level (N x N)
        SkylineMatrix	U;	    // U factor for the level (M x N)
        uvector<double> dval;	// reciprocal of the diagonal of the diagonal block

        uvector<uint>	map;	// indices map: permuted -> original
        uvector<uint>	rmap;	// indices map: original -> permuted

        // Misc variables
        mutable
        Vector r, w, x2, u0, F;	// some auxilary vectors for internal iterations

        Level() : q(0.0), niter(0), eps(0.0) { N = nnz = M = Md = 0; }
    };

    enum {
        INNER_ITER_FIXED,
        INNER_ITER_DYNAMIC
    } inner_iter_type;
    enum {
        LU_EXACT,
        LU_ILUT
    } lu_type;
    enum {
        ORDER_NONE,
        ORDER_ORIGINAL,
        ORDER_BLOCK,
        ORDER_SIMPLE_1
    } order;

    uint   ilut_p;
    double ilut_tau;
    uint   degree_max;
    double level_ratio;
    uint   coarse_n;

    uint nlevels;
    std::vector<Level> levels;

    const SkylineMatrix& level0_A;  // Matrix for level 0 is not kept in levels so that we don't
                                    // need to make a copy

#ifdef HAVE_UMFPACK
    // UMFPACK factorization of the coarsest level
    mutable void *Ac_symbolic, *Ac_numeric;
#endif

#ifdef PRINT_NORMS
    mutable std::ostringstream *norm_oss;
    void dump_norm_trace() const;
#endif

    void construct_level(uint level, const SkylineMatrix& A);

    void solve(uint level, const Vector& f, Vector& x) const;

    void solve_U(uint level, const Vector& w, Vector& x) const;
    void solve_L(uint level, const Vector& f, Vector& w) const;
    void solve_diagonal(uint level, const Vector& f, Vector& x) const;

    void order_none    (const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
                        const uvector<int>& nlinks_in, const uvector<int>& nlinks_out,
                        uint& Md, uint& M, uvector<uint>& map) const;
    void order_original(const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
                        uvector<int>& nlinks_in, uvector<int>& nlinks_out,
                        uint& Md, uint& M, uvector<uint>& map) const;
    void order_simple_1(const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
                        const uvector<int>& nlinks_in, const uvector<int>& nlinks_out,
                        uint& Md, uint& M, uvector<uint>& map) const;
    void order_block   (const SkylineMatrix& A, const LinkTypeMultiSplit& ltype_, const uvector<double>& aux,
                        const uvector<int>& nlinks_in, const uvector<int>& nlinks_out,
                        uint& Md, uint& M, uvector<uint>& map) const;

public:
    MultiSplitPrec(const SkylineMatrix& A, const Config& cfg);
    ~MultiSplitPrec();

    void solve(Vector& f, Vector& x) const; // Wrapper for solve(level,f,x)
    void graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const;

    friend std::ostream& operator<<(std::ostream& os, const MultiSplitPrec& p);
    friend std::ostream& operator<<(std::ostream& os, const MultiSplitPrec::Level& li);
};

#endif
