#include "matrix.h"
#include "include/logger.h"
#include "include/exception.h"
#include "include/tools.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <algorithm>

DEFINE_LOGGER("CSRMatrix");

CSRMatrix::CSRMatrix() {
    nrow = ncol = 0;
}

const CSRMatrix& CSRMatrix::operator=(const CSRMatrix& A) {
    nrow = A.nrow;
    ncol = A.ncol;
    mode = A.mode;

    if (nrow) {
	uint n = A.ia.size();
	uint nnz = A.ia[n-1];
	ia.resize(n);
	ja.resize(nnz);
	a.resize(nnz);
	memcpy(&ia[0], &A.ia[0], n*sizeof(int));
	memcpy(&ja[0], &A.ja[0], nnz*sizeof(int));
	memcpy(&a[0],  &A.a[0],  nnz*sizeof(double));

    } else {
	ia.clear();
	ja.clear();
	a.clear();
    }

    return *this;
}

double CSRMatrix::get(uint i, uint j) const THROW {
    // LOG_DEBUG("i = " << i << ", j = " << j);
    check_indices(i, j);

    std::vector<uint>::const_iterator start = ja.begin() + ia[i]; 
    std::vector<uint>::const_iterator   end = ja.begin() + ia[i+1];
    std::vector<uint>::const_iterator it = std::lower_bound(start, end, j);
    if (it != end && !(j < *it)) 
	return a[it-ja.begin()];

    LOG_WARN("Returning zero element for i = " << i << ", j = " << j);
    return 0;
}

void CSRMatrix::add(uint i, uint j, double v) THROW {
    check_indices(i, j);

    // in case matrix is in skyline form
    std::vector<uint>::iterator start = ja.begin() + ia[i]; 
    std::vector<uint>::iterator   end = ja.begin() + ia[i+1];
    std::vector<uint>::iterator it = std::lower_bound(start, end, j);
    if (it == end || j < *it) {
	// creating new element
	uint pos = it - ja.begin();
	ja.insert(it, j);
	a.insert (a.begin() + pos, v);
	for (uint k = i+1; k <= nrow; k++)
	    ia[k]++;
    } else {
	// adding to existing element
	a[it - ja.begin()] += v;
    }
}

void transpose(const CSRMatrix& A, CSRMatrix& B) {
    B.nrow = A.ncol;
    B.ncol = A.nrow;

    // if B has already contained matrix, flush it
    B.ia.clear();
    B.ja.clear();
    B.a.clear();

    std::vector<int> marker(A.ncol, 0);

    // Calculate number of elements in each column
    for (uint i = 0; i < A.nrow; i++)
	for (uint j = A.ia[i]; j < A.ia[i+1]; j++) 
	    marker[A.ja[j]]++;

    // Fill array ia for transposed matrix using marker
    // After cycle marker = B.ia (but shorter by 1 element)
    B.ia.resize(B.nrow+1);
    B.ia[0] = 0;
    for (uint i = 0; i < B.nrow; i++) {
	B.ia[i+1] = B.ia[i] + marker[i];
	marker[i] = B.ia[i];
    }
  
    // Fill in B.ja and B.a
    B.ja.resize(B.ia[B.nrow]);
    B.a.resize(B.ia[B.nrow]);
    int col, relpos;
    for (uint i = 0; i < A.nrow; i++)
	for (uint j = A.ia[i]; j < A.ia[i+1]; j++) {
	    col = A.ja[j];
	    relpos = marker[col]; 

	    B.ja[relpos] = i;
	    B.a[relpos] = A.a[j];

	    marker[col]++;
	}
}

bool CSRMatrix::is_symmetric() const {
    ASSERT(nrow == ncol, "Can not call is_symmetric() for non square matrices");
    for (uint i = 0; i < nrow; i++)
	for (uint j = ia[i]; j < ia[i+1]; j++)
	    if (a[j] != get(ja[j], i))
		return false;
    return true;
}

std::ostream& operator<<(std::ostream& os, const CSRMatrix& sm) {
    os << "Size: " << sm.nrow << "x" << sm.ncol << std::endl;
    for (uint i = 0; i < sm.nrow; i++) {
	os << "  Row: " << i << std::endl;
	for (uint j = sm.ia[i]; j < sm.ia[i+1]; j++)
	    os << "    " << sm.ja[j] << ": " << sm.a[j] << std::endl;
	os << std::endl;
    }

    return os;
}

void multiply(const CSRMatrix& A, const Vector& v, Vector& res, char type) THROW {
    ASSERT(res.size(), "Memory for result must have been already allocated");

    memset(&res[0], 0, res.size()*sizeof(double));
    switch (type) {
	case 'o': // A*v
	    ASSERT(A.ncol == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.nrow, "Not enough space in res vector");

	    for (uint i = 0; i < A.nrow; i++) 
		for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
		    res[i] += A.a[j] * v[A.ja[j]];
	    break;

	case 't': // A^T*v
	    ASSERT(A.nrow == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.ncol, "Not enough space in res vector");

	    for (uint i = 0; i < A.nrow; i++) 
		for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
		    res[A.ja[j]] += A.a[j] * v[i];
	    break;

	default:
	    THROW_EXCEPTION("Unknown type");
    }
}

void residual(const CSRMatrix& A, const Vector& b, const Vector& x, Vector& r) THROW {
    const int n = A.size();
    ASSERT((int)b.size() == n && (int)x.size() == n && (int)r.size() == n, "Wrong sizes");

    multiply(A, x, r);
    for (int i = 0; i < n; i++)
	r[i] = b[i] - r[i];
}
