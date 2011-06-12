#include "include/logger.h"
#include "include/tools.h"
#include "modules/solvers/solvers.h"
#include "bgs_prec.h"

#include <algorithm>

DEFINE_LOGGER("BGSPrec");

// #define lN (2*2)
#define lN (60*220)

void BGSPrec::solve(Vector& f, Vector& x) const THROW {
    ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    Vector f1(n1), f2(n2);
    for (uint z = 0; z*lN < n; z++)
	if (!(z&1)) memcpy(&f1[(z/2)*lN], &f[z*lN], lN*sizeof(double));
	else	    memcpy(&f2[(z/2)*lN], &f[z*lN], lN*sizeof(double));

    SolverStats stats;

    /* Solve A1 system */
    Vector x1(n1);
#ifndef PREC_SUBST
    DirectSolver(A1, f1, x1, A1_symbolic, A1_numeric, stats);
#else
    // SimpleSolver(A1, f1, *B1, x1, stats, 1e-2, NORM_L2, true);
    B1->solve(f1, x1);
#endif

    /* Solve A2 system */
    Vector r(n2);
    residual(A21, f2, x1, r);
    Vector x2(n2);
#ifndef PREC_SUBST
    DirectSolver(A2, r, x2, A2_symbolic, A2_numeric, stats);
#else
    // SimpleSolver(A2, r, *B2, x2, stats, 1e-2, NORM_L2, true);
    B2->solve(r, x2);
#endif

    for (uint z = 0; z*lN < n; z++)
	if (!(z&1)) memcpy(&x[z*lN], &x1[(z/2)*lN], lN*sizeof(double));
	else	    memcpy(&x[z*lN], &x2[(z/2)*lN], lN*sizeof(double));
}

BGSPrec::BGSPrec(const SkylineMatrix& A) {
    ASSERT(A.cols() == A.rows(), "Matrix must be square");
    n = A.rows();

    const uvector<uint>&  ia = A.get_ia();
    const uvector<uint>&  ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    uvector<uint>& a1_ia = A1.get_ia(), &a2_ia = A2.get_ia(), &a21_ia = A21.get_ia();
    uvector<uint>& a1_ja = A1.get_ja(), &a2_ja = A2.get_ja(), &a21_ja = A21.get_ja();
    uvector<double>& a1_a = A1.get_a(), &a2_a = A2.get_a(),   &a21_a = A21.get_a();

    a1_ia.push_back(0);
    a2_ia.push_back(0);
    a21_ia.push_back(0);

    a1_ja.resize(ja.size());
    a2_ja.resize(ja.size());
    a21_ja.resize(ja.size());
    a1_a.resize(a.size());
    a2_a.resize(a.size());
    a21_a.resize(a.size());

    uint ind1 = 0, ind2 = 0, ind21 = 0;
    uint i1 = 0, i2 = 0;
    for (uint i0 = 0, z = 0; i0 < n; i0 += lN, z++) {
	uint i = i0;

	if (!(z&1)) {
	    /* Block A1 */
	    uint i10 = i1;
	    for (uint k = 0; k < lN; k++, i1++, i++) {
		for (uint j = ia[i]; j < ia[i+1]; j++)
		    if (i0 <= ja[j] && ja[j] < i0+lN) {
			a1_ja[ind1] = i10 + (ja[j] - i0);
			a1_a[ind1] = a[j];
			ind1++;
		    }
		a1_ia.push_back(ind1);
	    }
	} else {
	    /* Blocks A2, A21 */
	    uint i20 = i2;
	    for (uint k = 0; k < lN; k++, i2++, i++) {
		for (uint j = ia[i]; j < ia[i+1]; j++)
		    if (i0 <= ja[j] && ja[j] < i0+lN) {
			a2_ja[ind2] = i20 + (ja[j] - i0);
			a2_a[ind2] = a[j];
			ind2++;
		    } else if (ja[j] < i0) {
			a21_ja[ind21] = i20 + (ja[j] - i0 + lN);
			a21_a[ind21] = a[j];
			ind21++;
		    }
		a2_ia.push_back(ind2);
		a21_ia.push_back(ind21);
	    }
	}
    }
    a1_ja.resize(a1_ia.back());
    a2_ja.resize(a2_ia.back());
    a21_ja.resize(a21_ia.back());
    a1_a.resize(a1_ia.back());
    a2_a.resize(a2_ia.back());
    a21_a.resize(a21_ia.back());

    A1.set_size(i1, i1);
    A2.set_size(i2, i2);
    A21.set_size(i2, i1);

    n1 = i1;
    n2 = i2;

    A1_symbolic = A1_numeric = NULL;
    A2_symbolic = A2_numeric = NULL;

#ifdef PREC_SUBST
    Config cfg;
    cfg.use_tails = true;
    cfg.niters = std_vector<uint>(2, 1, 1);
    cfg.sigmas.push_back(0.98);

#if 0
    B1 = new MultiSplitPrec(A1, cfg);
    B2 = new MultiSplitPrec(A2, cfg);
#else
    B1 = new DiagPrec(A1);
    B2 = new DiagPrec(A2);
#endif
#endif
}
