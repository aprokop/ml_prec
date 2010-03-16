#include "include/logger.h"
#include "include/tools.h"
#include "sym.h"
#include "modules/solvers/solvers.h"

DEFINE_LOGGER("SymPrec");

void SymPrec::solve(Vector& f, Vector& x) THROW {
    ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    PCGSolver(Asym, f, *B, x, eps, true);
}

SymPrec::SymPrec(const SkylineMatrix& A, const Config& cfg) {
    ASSERT(A.cols() == A.rows(), "Matrix must be square");

    n = A.size();

    Asym = A;
    const uvector<uint>& ia = Asym.get_ia();
    const uvector<uint>& ja = Asym.get_ja();
    uvector<double>&      a = Asym.a;

    /* Create a symmetric variant of the matrix */
    for (uint i = 0; i < n; i++)
	for (uint j = ia[i]+1; j < ia[i+1]; j++) {
	    double d = Asym(ja[j],i); /* <= 0 */
	    if (d > a[j])
		a[j] = d;
	}

    Config cfg1 = cfg;
    cfg1.unsym_matrix = false;

    B.reset(new Prec(Asym, cfg1));

    eps = 1e-14;
}
