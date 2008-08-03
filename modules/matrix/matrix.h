#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "config/config.h"
#include "include/exception.h"
#include "modules/vector/vector.h"

#include <map>

class MatrixInterface {
protected:
    uint nrow, ncol;

    void check_indices(uint i, uint j) const THROW {
	ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
	ASSERT(j < ncol, "Col index is out of boudaries: j = " << j << ", ncol = " << ncol);
    }

public:
    virtual ~MatrixInterface() {}
    virtual Vector operator*(const Vector& v) const THROW = 0;

    uint rows() const {
	return nrow;
    }
    uint cols() const {
	return ncol;
    }
    uint size() const THROW {
	ASSERT(nrow == ncol, "size was called for rectangular matrix: nrow = " << nrow << ", ncol = " << ncol);
	return nrow;
    }
};

// Class implements sparse matrix and operations with it
class CSRMatrix;
class SparseMatrix : public MatrixInterface {
protected:
    typedef std::map<uint, double> Row;
    std::vector<Row> vrows;

public:
    SparseMatrix(uint m = 0) THROW;
    SparseMatrix(uint m, uint n) THROW;
    SparseMatrix(const CSRMatrix& A) THROW;

    double  get(uint i, uint j) const THROW;
    void    set(uint i, uint j, double x) THROW;
    void    add(uint i, uint j, double x) THROW;

    Vector operator*(const Vector& v) const THROW {
	Vector x(nrow);
	multiply(*this, v, x);
	return x;
    }

    // Friend classes
    friend class CSRMatrix;

    // Friend functions
    friend void multiply(const SparseMatrix& A, const Vector& v, Vector& res, char type = 'o') THROW;
    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& sm);
};

class FVSparseMatrix : public SparseMatrix {
public:
    FVSparseMatrix(uint m = 0) THROW : SparseMatrix(m) { }
    FVSparseMatrix(uint m, uint n) THROW : SparseMatrix(m, n) { }

    void    new_link(uint i0, uint i1, double x) THROW;
    double  remove_link(uint i0, uint i1) THROW;
};

class CSRMatrix : public MatrixInterface {
protected:
    std::vector<uint> ia, ja;
    std::vector<double> a;

public:
    CSRMatrix();

    Vector operator*(const Vector& v) const THROW {
	Vector x(nrow);
	multiply(*this, v, x);
	return x;
    }
    double get(uint i, uint j) const THROW;

    bool is_symmetric() const;

    // Friend classes
    friend class SparseMatrix;
    friend class Mesh;
    friend class Prec;

    // Friend functions
    friend void	multiply(const CSRMatrix& A, const Vector&v, Vector& res, char type = 'o') THROW;
    friend void transpose(const CSRMatrix& A, CSRMatrix& B);
    friend std::ostream& operator<<(std::ostream& os, const CSRMatrix& sm);
};

typedef CSRMatrix SkylineMatrix;

void multiply(const MatrixInterface& A, const Vector& v, Vector& res, char type = 'o') THROW;

#endif // #ifndef __MATRIX_H__
