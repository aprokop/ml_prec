#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "config/config.h"
#include "include/exception.h"
#include "modules/vector/vector.h"

typedef unsigned int uint;

// We always assume that each row of CSRMatrix contains at least one element
class CSRMatrix {
protected:
    uint nrow, ncol;

    void check_indices(uint i, uint j) const THROW {
	ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
	ASSERT(j < ncol, "Col index is out of boudaries: j = " << j << ", ncol = " << ncol);
    }

    std::vector<uint> ia, ja;
    std::vector<double> a;

    enum {
	csr,
	skyline
    } mode;

public:
    CSRMatrix();
    virtual ~CSRMatrix() { }

    const CSRMatrix& operator=(const CSRMatrix& A);

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

    Vector operator*(const Vector& v) const THROW {
	Vector x(nrow);
	multiply(*this, v, x);
	return x;
    }
    virtual double get(uint i, uint j) const THROW;
    virtual void   add(uint i, uint j, double x) THROW;

    virtual bool is_symmetric() const;

    // Friend classes
    friend class SPEMesh;
    friend class Prec;
    friend class RelPrec;
    friend class AMGPrec;

    // Friend functions
    friend int main(int argc, char* argv[]);
    friend void	multiply(const CSRMatrix& A, const Vector& v, Vector& res, char type = 'o') THROW;
    friend void transpose(const CSRMatrix& A, CSRMatrix& B);
    friend std::ostream& operator<<(std::ostream& os, const CSRMatrix& sm);
};

class SkylineMatrix : public CSRMatrix {
public:
    SkylineMatrix();

    double get(uint i, uint j) const THROW;
    void   add(uint i, uint j, double x) THROW;

    bool exist(uint i, uint j) const THROW;
};

// r = b - Ax
void residual(const CSRMatrix& A, const Vector& b, const Vector& x, Vector& r) THROW;

#endif // #ifndef __MATRIX_H__
