#include "include/logger.h"
#include "include/tools.h"
#include "gs_prec.h"

#include <algorithm>

DEFINE_LOGGER("GSPrec");

void GSPrec::solve(Vector& f, Vector& x) const THROW {
    ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    const uvector<uint>&  ia = L.get_ia();
    const uvector<uint>&  ja = L.get_ja();
    const uvector<double>& a = L.get_a();

    for (uint i = 0; i < n; i++) {
	double d = f[i];
	for (uint j = ia[i]+1; j < ia[i+1]; j++)
	    d -= a[j]*x[ja[j]];

	x[i] = d/a[ia[i]];
    }

}

GSPrec::GSPrec(const SkylineMatrix& A) {
    ASSERT(A.cols() == A.rows(), "Matrix must be square");
    n = A.rows();

    const uvector<uint>&  ia = A.get_ia();
    const uvector<uint>&  ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    uvector<uint>&  lia = L.get_ia();
    uvector<uint>&  lja = L.get_ja();
    uvector<double>& la = L.get_a();

    lia.resize(ia.size());
    lja.resize(ja.size());
    la.resize(a.size());

    uint ind = 0;

    lia[0] = 0;
    for (uint i = 0; i < n; i++) {
	uint start = ia[i], end = ia[i+1];

	ASSERT(A(i,i) > 0, "Something is very strange: (" << i << ") has " << A(i,i) << " on diagonal, aborting...");

	/* Set diagonal element */
	lja[ind] = i;
	la[ind]  = a[start];
	ind++;

	uint middle = (std::upper_bound(ja.begin() + start+1, ja.begin() + end, i) - ja.begin());
	for (uint j = start+1; j < middle; j++) {
	    lja[ind] = ja[j];
	    la[ind]  = a[j];
	    ind++;
	}

	lia[i+1] = ind;
    }
    lja.resize(lia[n]);
    lja.resize(lia[n]);

    L.set_size(n,n);
}
