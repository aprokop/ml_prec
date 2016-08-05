#include "rel_prec.h"
#include "include/logger.h"

#include <iomanip>

DEFINE_LOGGER("Prec");

void RelPrec::solve(Vector& f, Vector& x) THROW {
    Vector _f(f);
    solve(_f, x, 0);
}

void RelPrec::solve(Vector& f, Vector& x, uint level) THROW {
    Level& li = levels[level];
    uint N = li.N;
    ASSERT(f.size() == N && x.size() == N, "Wrong dimension: N = " << N << ", f = " << f.size() << ", x = " << x.size());

    const std::vector<uint>& dtr   = li.dtr;
    const std::vector<uint>& tr    = li.tr;
    const std::vector<Tail>& tails = li.tails;

    // truncate tails
    for (uint i = 0; i < tails.size(); i++) {
        const Tail& tail = tails[i];
        f[tail[0].index] *= tail[0].a2;

        for (uint j = 1; j < tail.size(); j++) {
            const TailNode& tn = tail[j];
            f[tn.index] = tn.a2*f[tn.index] + tn.a3*f[tail[j-1].index];
        }
    }

    if (level < nlevels-1) {
        uint n = levels[level+1].N;

        Vector& f1 = li.f1;
        for (uint i = 0; i < n; i++)
            f1[i] = f[tr[i]];

        Vector& x1 = li.x1;
        if (niter) {
            const CSRMatrix& A = levels[level+1].A;

            Vector& tmp0 = li.tmp0;
            Vector& tmp1 = li.tmp1;

            tmp0 = f1;
            // ===============    STEP 1    ===============
            // x1 = B^{-1}f
            solve(tmp0, x1, level+1);

            // ===============    STEP 2+    ===============
            for (uint i = 1; i < niter; i++) {
                // x1 = x0 + B^{-1}(f - A*x0)
                residual(A, f1, x1, tmp0);
                solve(tmp0, tmp1, level+1);
                daxpy(1., tmp1, x1);
            }
        } else {
            solve(f1, x1, level+1);
        }

        for (uint i = 0; i < n; i++)
            x[tr[i]] = x1[i];

    } else {
        // for last level assert for now that we have only diagonal matrix
        ASSERT(tr.size() == 0, "We must have an \"almost\" diagonal on the coarsest level");
    }

    // solve diagonal part
    for (uint i = 0; i < dtr.size(); i++) {
        uint ii = dtr[i];
        x[ii] = f[ii] / li.aux[ii];
    }

    // restore tails in *backward* order
    for (int i = tails.size()-1; i >= 0; i--) {
        const Tail& tail = tails[i];
        if (tail.end_type == 'f')
            x[tail.back().index] = f[tail.back().index];
        for (int j = tail.size() - 2; j >= 0; j--) {
            const TailNode& tn = tail[j];
            x[tn.index] = tn.a1*x[tail[j+1].index] + f[tn.index];
        }
    }
}
