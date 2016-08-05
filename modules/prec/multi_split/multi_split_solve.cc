#include "include/logger.h"
#include "multi_split_prec.h"
#include "project/config.h"

#include "modules/solvers/solvers.h"

#include <iomanip>

DEFINE_LOGGER("MultiSplitPrec");

// Solve diagonal subsystem
void MultiSplitPrec::solve_diagonal(uint level, const Vector& f, Vector& x) const THROW {
    const Level& li             = levels[level];
    const uvector<uint>& map	= li.map;
    const uvector<double>& dval = li.dval;
    uint Md = li.Md;
    uint N  = li.N;
    for (uint i = 0; i < Md; i++) {
        uint j = map[(N-Md) + i];
        x[j] = dval[i]*f[j];
    }
}

// Solve L*w = f
void MultiSplitPrec::solve_L(uint level, const Vector& f, Vector& w) const {
    const Level& li = levels[level];
    uint N = li.N, Md = li.Md;
    const uvector<uint>& map = li.map;
    const CSRMatrix& L = li.L;

    for (uint i = 0; i < N-Md; i++) { /* i is a permuted index */
        w[i] = f[map[i]];
        for (uint j_ = L.ia[i]; j_ < L.ia[i+1]; j_++)
            w[i] -= L.a[j_] * w[L.ja[j_]];
    }
}

// Solve U*x = w
void MultiSplitPrec::solve_U(uint level, const Vector& w, Vector& x) const {
    const Level& li = levels[level];
    uint N = li.N, Md = li.Md, M = li.M;
    const uvector<uint>& map = li.map;
    const SkylineMatrix& U = li.U;

    for (int i = M-1; i >= 0; i--) {
        double z = w[i];

        for (uint j_ = U.ia[i]+1; j_ < U.ia[i+1]; j_++)
            z -= U.a[j_] * x[map[U.ja[j_]]];

        x[map[i]] = z / U.a[U.ia[i]];
    }
}

void MultiSplitPrec::solve(Vector& f, Vector& x) const THROW {
    solve(0, f, x);
}

void MultiSplitPrec::solve(uint level, const Vector& f, Vector& x) const THROW {
    const Level& li = levels[level];

#ifdef PRINT_NORMS
    /* Log norm */
    if (level < 3)
        (*norm_oss) << level << " " << dnrm2(f) << std::endl;
#endif

    uint N = li.N;
    uint M = li.M;
    uint n = (level < nlevels-1) ? levels[level+1].N : 0;

    uint niter = li.niter;

    const CSRMatrix& A = level ? levels[level].A : level0_A;
    const uvector<uint>& map = li.map;

    if (level == nlevels-1 && coarse_n) {
        SolverStats stats;
        DirectSolver(A, f, x, Ac_symbolic, Ac_numeric, stats);
        return;
    }

    Vector& r  = li.r;
    Vector& w  = li.w;
    Vector& u0 = li.u0;
    Vector& x2 = li.x2;
    Vector& F  = li.F;

    /* ===============    STEP 1    =============== */
    /*
     * x1 = x0 + solve(f1 - A*x0)
     * But x0 = 0, so x1 = solve(f1)
     */
    r = f;
    double init_norm = 0;
    if (li.eps)
        init_norm = dnrm2(r);
    solve_L(level, r, w);
    if (level < nlevels-1) {
        memcpy(&F[0], &w[M], n*sizeof(double));

        solve(level+1, F, x2);

        for (uint i = 0; i < n; i++)
            x[map[i+M]] = x2[i];
    }
    solve_U(level, w, x);
    solve_diagonal(level, r, x);

    if (level == nlevels-1) {
        residual(A, f, x, r);
        LOG_DEBUG("cnorm = " << dnrm2(r));
    }

    /* ===============    STEP 2+    =============== */
    /* x^{k+1} = x^k + solve(f1 - A*x^k) */
    for (uint i = 2; li.eps || i <= niter; i++) {
        residual(A, f, x, r);

#ifdef PRINT_NORMS
        /* Log norm */
        if (level < 3)
            (*norm_oss) << level << " " << dnrm2(r) << std::endl;
#endif

        if (li.eps) {
            double norm = dnrm2(r);
            if (norm / init_norm < li.eps) {
                LOG_INFO("Level #" << level << ": Number of level iterations: " << i);
                break;
            }
        }

        solve_L(level, r, w);
        if (level < nlevels-1) {
            memcpy(&F[0], &w[M], n*sizeof(double));

            solve(level+1, F, x2);

            for (uint i = 0; i < n; i++)
                u0[map[i+M]] = x2[i];
        }
        solve_U(level, w, u0);
        solve_diagonal(level, r, u0);

        daxpy(1., u0, x);
    }
#ifdef PRINT_NORMS
    /* Log norm */
    if (0 < level && level < 3) {
        residual(A, f, x, r);
        (*norm_oss) << level << " " << dnrm2(r) << std::endl;
    }
#endif
}
