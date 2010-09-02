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
    double v;
    for (uint i = 0; i < n; i++) {
	v = A(i,i);
	if (v <= 0)
	    THROW_EXCEPTION("Something is very strange: (" << i << ") has " << v << " on diagonal, aborting...");
	d[i] = 1./v;
    }
}
