#ifndef __CHEB_PREC_H__
#define __CHEB_PREC_H__

#include "config/config.h"
#include "include/time.h"
#include "modules/prec/misc/misc.h"
#include "project/config.h"
#include "modules/prec/prec_base.h"
#include "modules/mesh/mesh.h"
#include <iostream>

class LinkTypeCheb;
class Prec : public PrecBase {
private:
    bool use_tails;

    struct Level {
        uint N, nnz;		        // number of nodes and nonzeros elements for the level
        uint M;                     // size of the excluded block (not including diagonal block)
        uint Md;		            // size of the excluded diagonal block

        double alpha, beta;	        // spectral constants used for the level
        double lmin, lmax;	        // spectral boundaries for the level
        uint ncheb;		            // number of Chebyshev iterations on the level

        SkylineMatrix A;	        // level matrix (for level 0 we use level0_A)
        CSRMatrix     L;	        // L factor for the level (N x N)
        SkylineMatrix U;	        // U factor for the level (M x N)

        uvector<uint> map;	        // indices map: permuted -> original
        uvector<uint> rmap;	        // indices map: original -> permuted

        uvector<double> dval;	    // reciprocal of the diagonal of the diagonal block

        mutable
        Vector w, tmp, x2, u1, u0;  // Some auxilary vectors for Chebyshev iterations and Schur

        uvector<double> aux;	    // Auxilary array, corresponds to value c of the node; also
                                    // used in tail elimination. Don't actually need to keep it
    };

    uint nlevels;
    std::vector<Level> levels;

    const SkylineMatrix& level0_A;  // Matrix for level 0 is not kept in levels so that we don't
                                    // need to make a copy

    // This function tries to increase s by amount necessary to mark v as
    // removable. If it can not be done (i.e. there is not enough value of c
    // for that then s is not changed
    //
    // NOTE: if matrix contains positive offdiagonal elements (as in the case
    // when we take Schur complement to u and p in mixed hybrid FEM) removing
    // negative offdiagonal elements may result in matrix which is not positive
    // definite. Sometimes we can assume that if c is big enough new matrix
    // would still be positive definite
    bool to_remove(double c, double v, double beta, double& s) {
#if 1
        /* Dynamic variant */
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
#else
        /* Static variant */
        if (1 + 12*(-v)/c > beta)
            return false;
#endif

        return true;
    }

    void construct_permutation(const SkylineMatrix& A, const LinkTypeCheb& ltype, uvector<int>& nlinks,
                               uint& Md, uint& M, uvector<uint>& map) const;

    /* ViTE framework */
    mutable std::ostringstream *foss;
    double vite_start;
#ifdef USE_VITE
    void log_state(const char* name) const {
        *foss << "10 " << pclock() - vite_start << " S T0 " << name << std::endl;
    }
    void log_event(uint val) const {
        *foss << "11 " << pclock() - vite_start << " E T0 " << val << std::endl;
    }
#else
    void log_state(const char* name) const { }
    void log_event(uint val) const { }
#endif
    void dump_vite_trace() const;

    mutable std::ostringstream *norm_oss;
    void dump_norm_trace() const;

    void optimize_storage();
    void construct_level(uint level, const SkylineMatrix& A);

    void solve(uint level, const Vector& f, Vector& x) const;

public:
    Prec(const SkylineMatrix& A, const Config& cfg);
    ~Prec();

    void solve(Vector& f, Vector& x) const; /* Wrapper for solve(level,f,x) */
    void graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const;

    /* If the outer procedure is Chebyshev, it needs these values */
    double lmin() const { return levels[0].lmin; }
    double lmax() const { return levels[0].lmax; }

    friend std::ostream& operator<<(std::ostream& os, const Prec& p);
    friend std::ostream& operator<<(std::ostream& os, const Prec::Level& li);
};

#endif // __CHEB_PREC_H__
