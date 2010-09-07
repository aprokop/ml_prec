#include "misc.h"
#include "include/uvector.h"
#include "include/svector.h"

/* Saad. Iterative methods for sparse linear systems. Pages 310-312 */
void construct_sparse_lu(const SkylineMatrix& A, const uvector<uint>& map, const uvector<uint>& rmap,
			 uint Md, uint M, const LinkTypeBase& ltype, double beta, const uvector<double>& aux,
			 SkylineMatrix& nA, SkylineMatrix& U, CSRMatrix& L) {
    const uint	N = A.size();
    uint	n = N-M-Md;

    uvector<int> jr(N, -1);
    typedef svector<uint> container;

    container jw;
    uvector<double> w;
    uint max_num;

    const uvector<uint>  &A_ia = A.get_ia();
    const uvector<uint>  &A_ja = A.get_ja();
    const uvector<double> &A_a = A.get_a();
    uvector<uint>  &L_ia = L.get_ia(), &U_ia = U.get_ia(), &nA_ia = nA.get_ia();
    uvector<uint>  &L_ja = L.get_ja(), &U_ja = U.get_ja(), &nA_ja = nA.get_ja();
    uvector<double> &L_a = L.get_a(),   &U_a = U.get_a(),   &nA_a = nA.get_a();

    L.set_size  (N-Md, N-Md);
    U.set_size  (   M, N-Md);
    nA.set_size(   n,    n);

    L_ia.push_back(0);
    U_ia.push_back(0);
    nA_ia.push_back(0);

    /* Reserve some space */
    L_ja.reserve(2*(N-Md));       L_a.reserve(2*(N-Md));
    U_ja.reserve(3*M);		  U_a.reserve(3*M);
    nA_ja.reserve(5*(N-M-Md));    nA_a.reserve(5*(N-M-Md));
    /* TODO: deal with M = 0 */
    for (uint i = 0; i < N-Md; i++) { /* i corresponds to a permuted index */
	/* Step 0: clear tmp values filled on previous iteration */
	unsigned jwn = jw.size();
	for (uint k = 0; k < jwn; k++)
	    jr[jw[k]] = -1;
	jw.clear();
	w.clear();

	uint arow = map[i];	/* Index of a row in A */

	/* Step 1: create buffer with permuted row of A */
	/* Add diagonal element to buffer */
	jr[i] = 0;
	jw.insert(i);
	w.push_back(aux[map[i]]);    /* w[0] is the value of the diagonal element */

	/* Add off-diagonal elements to buffer */
	max_num = 1;
	ltype.set_row(arow);
	for (uint j_ = A_ia[arow]+1; j_ < A_ia[arow+1]; j_++) {
	    uint j = A_ja[j_];

	    if (ltype.stat(j_) == PRESENT) {
		/* Translate original indices into permuted */
		uint new_j = rmap[j];   /* permuted index */

		jr[new_j] = max_num++;
		jw.insert(new_j);

		/* Scale elements: required by our preconditioner theory. See report */
		double z = A_a[j_] / beta;

		w.push_back(z);
		w[0] -= z;       /* update diagonal */
	    }
	}

	/* Step 2: perform sparse gaussian elimination */
	uint m = std::min(i, M);

	uint k = *(jw.begin());
	while (k < m) {
	    double lik = w[jr[k]] / U(k,k);

	    /* Update L */
	    L_ja.push_back(k);
	    L_a.push_back(lik);

	    /* Update buffer using k-th row of U */
	    for (uint j_ = U_ia[k]+1; j_ < U_ia[k+1]; j_++) {
		uint j = U_ja[j_];      /* j is a permuted index */

		if (jr[j] != -1) {
		    /* Element already exists in the buffer => update */
		    w[jr[j]] -= lik*U_a[j_];
		} else {
		    /* Element does not exist in the buffer => create */
		    jr[j] = max_num++;
		    jw.insert(j);
		    w.push_back(-lik*U_a[j_]);
		}
	    }

	    /* Find next index in the buffer, i.e. next element of L */
	    k = *(jw.upper_bound(k));
	}
	L_ia.push_back(L_ja.size());

	/* Step 3: move element from buffer to corresponding rows of U/A_{level+1} */
	if (i < M) {
	    /* Update U */
	    for (container::const_iterator it = jw.lower_bound(i); it != jw.end(); it++) {
		U_ja.push_back(*it);
		U_a.push_back(w[jr[*it]]);
	    }

	    U_ia.push_back(U_ja.size());
	} else {
	    /* Update A_{level+1}
	     * Note that the process is a bit more difficult than updating U as A is a SkylineMatrix
	     * but elements added are sorted (i.e. the diagonal element is somewhere in the middle */
	    uint adind = nA_ja.size();
	    nA_ja.push_back(i-M);
	    nA_a.push_back(w[0]);

	    for (container::const_iterator it = jw.lower_bound(M); it != jw.end(); it++)
		if (*it != i) {
		    nA_ja.push_back(*it-M);
		    nA_a.push_back(w[jr[*it]]);
		}

	    nA_ia.push_back(nA_ja.size());
	}
    }
}
