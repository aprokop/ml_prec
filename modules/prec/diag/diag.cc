#include "include/logger.h"
#include "include/tools.h"
#include "diag_prec.h"

DEFINE_LOGGER("DiagPrec");

void DiagPrec::solve(Vector& f, Vector& x) const THROW {
    ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    for (uint i = 0; i < n; i++)
        x[i] = f[i] * d[i];
}

DiagPrec::DiagPrec(const SkylineMatrix& A) {
    ASSERT(A.cols() == A.rows(), "Matrix must be square");
    n = A.rows();

    d.resize(n);
    for (uint i = 0; i < n; i++) {
        ASSERT(A(i,i) > 0, "Something is very strange: (" << i << ") has " << A(i,i) << " on diagonal, aborting...");
        d[i] = 1./A(i,i);
    }
}
