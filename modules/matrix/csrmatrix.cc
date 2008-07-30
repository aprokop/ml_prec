#include "matrix.h"
#include "include/logger.h"
#include "include/tools.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

DEFINE_LOGGER("CSRMatrix");

CSRMatrix::CSRMatrix() {
    nrow = ncol = 0;
}

CSRMatrix::CSRMatrix(const SparseMatrix& A) {
    *this = A;
}

void CSRMatrix::operator=(const SparseMatrix& A) {
    nrow = A.rows();
    ncol = A.cols();
    
    // reserve for 2D
    ia.resize(nrow + 1, 0);
    ja.reserve(nrow*5);
    a.reserve(nrow*5);

    int ind = 0;
    for (int i = 0; i < nrow; i++) {
	const SparseMatrix::Row& row = A.vrows[i];
	ia[i+1] = ia[i] + row.size();
	for (SparseMatrix::Row::const_iterator it = row.begin(); it != row.end(); it++) {
	    ja.push_back(it->first);
	    a.push_back(it->second);
	    ind++;
	}
    }
}

void multiply(const CSRMatrix& A, const Vector&v, Vector& res, char type) THROW {
    ASSERT(res.size(), "Memory for result must have been already allocated");

    memset(&res[0], 0, res.size()*sizeof(double));
    switch (type) {
	case 'o': // A*v
	    ASSERT(A.ncol == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.nrow, "Not enough space in res vector");

	    for (int i = 0; i < A.nrow; i++) 
		for (int j = A.ia[i]; j < A.ia[i+1]; j++)
		    res[i] += A.a[j] * v[A.ja[j]];
	    break;

	case 't': // A^T*v
	    ASSERT(A.nrow == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.ncol, "Not enough space in res vector");

	    for (int i = 0; i < A.nrow; i++) 
		for (int j = A.ia[i]; j < A.ia[i+1]; j++)
		    res[A.ja[j]] += A.a[j] * v[i];
	    break;

	default:
	    THROW_EXCEPTION("Unknown type");
    }
}

double CSRMatrix::get(uint i, uint j) const THROW {
    check_indices(i, j);

    std::vector<int>::const_iterator start = ja.begin() + ia[i], end = ja.begin() + ia[i+1];
    std::vector<int>::const_iterator it = std::lower_bound(start, end, j);
    if (it != end && !(*it < j)) 
	return a[it-ja.begin()];

    LOG_WARN("Returning zero element for i = " << i << ", j = " << j);
    return 0;
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
    for (int i = 0; i < A.nrow; i++)
	for (int j = A.ia[i]; j < A.ia[i+1]; j++) 
	    marker[A.ja[j]]++;

    // Fill array ia for transposed matrix using marker
    // After cycle marker = B.ia (but shorter by 1 element)
    B.ia.resize(B.nrow+1);
    B.ia[0] = 0;
    for (int i = 0; i < B.nrow; i++) {
	B.ia[i+1] = B.ia[i] + marker[i];
	marker[i] = B.ia[i];
    }
  
    // Fill in B.ja and B.a
    B.ja.resize(B.ia[B.nrow]);
    B.a.resize(B.ia[B.nrow]);
    int col, relpos;
    for (int i = 0; i < A.nrow; i++)
	for (int j = A.ia[i]; j < A.ia[i+1]; j++) {
	    col = A.ja[j];
	    relpos = marker[col]; 

	    B.ja[relpos] = i;
	    B.a[relpos] = A.a[j];

	    marker[col]++;
	}
}

std::ostream& operator<<(std::ostream& os, const CSRMatrix& sm) {
    os << "Size: " << sm.nrow << "x" << sm.ncol << std::endl;
    for (int i = 0; i < sm.nrow; i++) {
	os << "  Row: " << i << std::endl;
	for (int j = sm.ia[i]; j < sm.ia[i+1]; j++)
	    os << "    " << sm.ja[j] << ": " << sm.a[j] << std::endl;
	os << std::endl;
    }

    return os;
}
