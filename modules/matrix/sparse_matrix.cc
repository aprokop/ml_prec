#include "matrix.h"
#include "include/logger.h"

#include <iostream>

DEFINE_LOGGER("SparseMatrix");

SparseMatrix::SparseMatrix(uint m) THROW {
    vrows.resize(m);
    nrow = ncol = m;
}

SparseMatrix::SparseMatrix(uint m, uint n) THROW {
    vrows.resize(m);
    nrow = m;
    ncol = n;
}

SparseMatrix::SparseMatrix(const CSRMatrix& A) THROW {
    nrow = A.rows();
    ncol = A.cols();

    vrows.resize(nrow);
    for (uint i = 0; i < nrow; i++) 
	for (uint j = A.ia[i]; j < A.ia[i+1]-1; j++) 
	    vrows[i][A.ja[j]] = A.a[j];
}

double SparseMatrix::get(uint i, uint j) const THROW {
    check_indices(i, j);

    Row::const_iterator it = vrows[i].find(j);
    if (it != vrows[i].end())
	return it->second;

    LOG_WARN("Returning zero element");
    return 0;
}

void SparseMatrix::set(uint i, uint j, double x) THROW {
    check_indices(i, j);

    vrows[i][j] = x;
}

void SparseMatrix::add(uint i, uint j, double x) THROW {
    check_indices(i, j);

    Row::iterator it = vrows[i].find(j);
    if (it != vrows[i].end())
	it->second += x;
    else
	vrows[i][j] = x;
}

void multiply(const SparseMatrix& A, const Vector&v, Vector& res, char type) THROW {
    switch (type) {
	case 'o': // A*v
	    ASSERT(A.ncol == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.nrow, "Not enough space in res vector");
	    memset(&res[0], 0, res.size()*sizeof(double));

	    for (uint i = 0; i < A.nrow; i++) 
		for (SparseMatrix::Row::const_iterator it = A.vrows[i].begin(); it != A.vrows[i].end(); it++)
		    res[i] += it->second * v[it->first];
	    break;

	case 't': // A^T*v
	    ASSERT(A.nrow == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.ncol, "Not enough space in res vector");
	    memset(&res[0], 0, res.size()*sizeof(double));

	    for (uint i = 0; i < A.nrow; i++) 
		for (SparseMatrix::Row::const_iterator it = A.vrows[i].begin(); it != A.vrows[i].end(); it++)
		    res[it->first] += it->second * v[i];
	    break;

	default:
	    THROW_EXCEPTION("Unknown type");
    }
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& sm) {
    os << "Size: " << sm.nrow << "x" << sm.ncol << std::endl;

    for (uint i = 0; i < sm.nrow; i++, os << std::endl) {
	os << "   ";
	for (uint j = 0; j < sm.ncol; j++)
	    os << sm.get(i,j) << " ";
    }

    return os;
}

void multiply(const MatrixInterface& A, const Vector& v, Vector& res, char type) THROW {
    try {
	multiply(dynamic_cast<const SparseMatrix&>(A), v, res, type);
    } catch (std::bad_cast) {
	multiply(dynamic_cast<const CSRMatrix&>(A), v, res, type);
    }
}
