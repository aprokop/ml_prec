#include "include/logger.h"
#include "include/tools.h"
#include "modules/solvers/solvers.h"
#include "rbgs_prec.h"

#include <algorithm>
#include <numeric>

DEFINE_LOGGER("RBGSPrec");

void RBGSPrec::solve(Vector& f, Vector& x) const THROW {
    ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    SolverStats stats;
#ifndef PREC_SUBST
    DirectSolver(B, f, x, B_symbolic, B_numeric, stats);
#else
    Bprec->solve(f, x);
#endif
}

RBGSPrec::RBGSPrec(const SkylineMatrix& A, const Config& cfg) {
    ASSERT(A.cols() == A.rows(), "Matrix must be square");
    n = A.rows();
    lN = cfg.nx*cfg.ny;

    const uvector<uint>&  ia = A.get_ia();
    const uvector<uint>&  ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    uvector<uint>   &b_ia = B.get_ia(), &b_ja = B.get_ja();
    uvector<double> &b_a = B.get_a();

    b_ia.resize(n+1);
    b_ja.resize(ja.size());
    b_a.resize(a.size());

    uint R_nnz = 0, R_nnz_total = 0;

    uint ind = 0;
    b_ia[0] = 0;
    for (uint i = 0; i < n; i++) {
        double c = std::accumulate(&a[ia[i]], &a[ia[i+1]], 0.0);
        for (uint j = ia[i]; j < ia[i+1]; j++)
            if (ja[j] <= i || ((ja[j] - i) % lN)) {
                b_ja[ind] = ja[j];
                b_a[ind]  = a[j];
                ind++;
            } else {
                double v = -a[j];
                if (v > 0.5*alpha*c/(1-alpha)) {
                    b_ja[ind] = ja[j];
                    b_a[ind]  = a[j];
                    ind++;

                    R_nnz++;
                }
                R_nnz_total++;
            }
        b_ia[i+1] = ind;
    }
    b_ja.resize(b_ia.back());
    b_a.resize(b_ia.back());

    LLL_INFO("R_nnz = " << R_nnz);
    LLL_INFO("R_nnz_total = " << R_nnz_total);

    B.set_size(n, n);

    B_symbolic = B_numeric = NULL;

#ifdef PREC_SUBST
    Config cfg;
    cfg.use_tails = true;
    cfg.niters = std_vector<uint>(2, 1, 5);
    cfg.sigmas.push_back(0.3);

#if 1
    Bprec = new MultiSplitPrec(B, cfg);
#else
    Bprec = new DiagPrec(B);
#endif
#endif
}
