#include "matrix.h"
#include "include/logger.h"
#include "include/exception.h"
#include "include/tools.h"

#include <algorithm>

DEFINE_LOGGER("SkylineMatrix");

SkylineMatrix::SkylineMatrix() {
    nrow = ncol = 0;
}

uint SkylineMatrix::index(uint i, uint j) const {
    check_indices(i, j);

    if (i == j)
	return ia[i];

    uvector<uint>::const_iterator start = ja.begin() + ia[i] + 1; 
    uvector<uint>::const_iterator   end = ja.begin() + ia[i+1];
    uvector<uint>::const_iterator it = std::lower_bound(start, end, j);
    if (it != end && !(j < *it)) 
	return it-ja.begin();

    return uint(-1);
}

void SkylineMatrix::load(const std::string& filename, bool transform, bool ascii) THROW {
    CSRMatrix::load(filename, ascii);

    if (transform) {
	/* Transform matrix from CSR to Skyline */
	std::vector<uint>::const_iterator it;
	uint start, end, dind = 0xffffffff, j;
	for (uint i = 0; i < nrow; i++) {
	    start = ia[i];
	    end   = ia[i+1];

	    for (j = start; j < end; j++)
		if (ja[j] == i) {
		    dind = j;
		    break;
		}
	    if (j == end)
		THROW_EXCEPTION("No diagonal element in row " << i);

	    double v = a[dind];
	    for (uint j = dind; j > start; j--) {
		ja[j] = ja[j-1];
		a[j]  = a[j-1];
	    }
	    ja[start] = i;
	    a[start] = v;
	}
    }
}

SkylineMatrix::SkylineMatrix(const MapMatrix& A) {
    nrow = A.rows();
    ncol = A.cols();
    ia.resize(nrow + 1);
    ja.reserve(nrow*7);
    a.reserve(nrow*7);

    ia[0] = 0;
    for (uint i = 0; i < nrow; i++) {
	const MapMatrix::Row& row = A.data[i];

	/* First we add diagonal element */
	ja.push_back(i);
	a.push_back(row.find(i)->second);

	/* Add remaining elements of a row */
	for (MapMatrix::Row::const_iterator it = row.begin(); it != row.end(); it++) 
	    if (it->first != i) {
		ja.push_back(it->first);
		a.push_back(it->second);
	    }

	ia[i+1] = ia[i] + row.size();
    }
}

void SkylineMatrix::optimize_storage(char type) {
    if (type == 's') {
	/* Add symmetric storage */
	uint n = size(), nnz = ja.size();
	uvector<uint>::const_iterator jastart = ja.begin();

	iasym.resize(n+1);
	jasym.resize((nnz+n) >> 1);
	asym.resize((nnz+n) >> 1);

	iasym[0] = 0;
	uint ind = 0;
	for (uint i = 0; i < nrow; i++) {
	    uint start = ia[i], end = ia[i+1];
	    uint middle = std::upper_bound(jastart + start+1, jastart + end, i) - jastart;

	    /* Store lower half */
	    for (uint _j = start; _j < middle; _j++, ind++) {
		jasym[ind] = ja[_j];
		asym[ind]  =  a[_j];
	    }

	    iasym[i+1] = ind;
	}
    } else {
	THROW_EXCEPTION("Unknown type: " << type);
    }
}

void sym_multiply(const SkylineMatrix& A, const Vector& v, Vector& res) THROW {
    ASSERT(A.rows() == res.size(), "Different sizes: A is " << A.sizes() << ", res is " << res.size());
    ASSERT(A.ncol == v.size(), "Multiplying sparse matrix and vector with different dimensions");
    ASSERT(res.size() == A.nrow, "Not enough space in res vector");

    memset(&res[0], 0, res.size()*sizeof(double));
    if (A.iasym.size()) {
	/* Matrices were optimized */
	for (uint i = 0; i < A.nrow; i++) {
	    uint start = A.iasym[i], end = A.iasym[i+1];

	    /* Add diagonal */
	    res[i] += A.asym[start]*v[i];

	    /* Rest of elements */
	    for (uint _j = start+1; _j < end; _j++) {
		uint j = A.jasym[_j];
		res[i] += A.asym[_j] * v[j];
		res[j] += A.asym[_j] * v[i];
	    }
	}
    } else {
	/* Matrices were NOT optimized */
	for (uint i = 0; i < A.nrow; i++) {
	    uint start = A.ia[i], end = A.ia[i+1];
	    uint middle = (std::upper_bound(A.ja.begin() + start+1, A.ja.begin() + end, i) - A.ja.begin());

	    /* Add diagonal */
	    res[i] += A.a[start]*v[i];

	    /* Rest of elements */
	    for (uint _j = start+1; _j < middle; _j++) {
		uint     j = A.ja[_j];
		double val = A.a[_j];

		res[i] += val * v[j];
		res[j] += val * v[i];
	    }
	}
    }
}

static void psort(const uint *a, uint n, std::vector<uint>& sorted) {
    int j;
    uint c;
    uint v;

    sorted[0] = 0;
    for (uint i = 1; i < n; i++) {
	c = sorted[i] = i;

	v = a[c];
	for (j = i-1; j >= 0 && a[sorted[j]] > v; j--)
	    sorted[j+1] = sorted[j];
	sorted[j+1] = c;
    }
}

void SkylineMatrix::permute(const std::vector<uint>& perm) THROW {
    uint n = size();
    ASSERT(perm.size() == n, "n = " << n << ", but perm.size() = " << perm.size());

    /* Compute inverse permutation vector */
    std::vector<uint> iperm(n);
    for (uint i = 0; i < n; i++) 
	iperm[perm[i]] = i;

    SkylineMatrix B = *this;

    /* Construct new ia */
    B.ia[0] = 0;
    for (uint bi = 0; bi < n; bi++) 
	B.ia[bi+1] = B.ia[bi] + (ia[iperm[bi]+1] - ia[iperm[bi]]);

    std::vector<uint> inds, sorted;
    for (uint bi = 0; bi < n; bi++) {
	uint ai = iperm[bi];
	uint rstart = B.ia[bi], rend = B.ia[bi+1];

	B.ja[rstart] = bi;
	B.a[rstart] = a[ia[ai]];

	uint rn = ia[ai+1] - ia[ai] - 1;

	inds.resize(rn);
	sorted.resize(rn);
	for (uint j = 0; j < rn; j++)
	    inds[j] = perm[ja[ia[ai] + j+1]];
	psort(&inds[0], rn, sorted);

	for (uint j = 0; j < rn; j++) {
	    B.ja[rstart+j+1] = inds[sorted[j]];
	    B.a[rstart+j+1]  = a[ia[ai]+1 + sorted[j]];
	}
    }

    (*this) = B;
}
