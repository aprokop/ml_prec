#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "config/config.h"
#include "include/exception.h"
#include "modules/vector/vector.h"
#include "include/uvector.h"

#include <map>
#include <cstring>

class Matrix {
protected:
    uint nrow;
    uint ncol;

    void check_indices(uint i, uint j) const THROW {
	ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
	ASSERT(j < ncol, "Col index is out of boudaries: j = " << j << ", ncol = " << ncol);
    }

public:
    virtual ~Matrix() {}
    // virtual Vector operator*(const Vector& v) const THROW = 0;

    uint rows() const { return nrow; }
    uint cols() const { return ncol; }

    virtual double  operator()(uint i, uint j) const THROW = 0;
    virtual double& operator()(uint i, uint j) THROW = 0;

    std::pair<uint,uint> sizes() const {
	return std::pair<uint,uint>(nrow,ncol);
    }
};

class CSRMatrix;
/* Dense matrix class */
class DMatrix : public Matrix {
private:
    std::vector<double> data;

    /*
     * Symmetric matrix can also be stored in factored form.
     * In this case the original matrix is stored in upper triangle,
     * L (except the diagonal) is stored in lower triangle, diagonal
     * of L is stored in d
     */
    bool factored;
    std::vector<double> d;

public:
    enum Type {
	ZERO,
	IDENTITY
    };

public:
    DMatrix(uint _n = 1, Type type = ZERO) THROW;
    DMatrix(uint _m, uint _n) THROW;
    DMatrix(const DMatrix& v);
    DMatrix(const char* values) THROW;
    DMatrix(const std::vector<double>& v, uint _ncol) THROW;
    DMatrix(const CSRMatrix& sm);
    ~DMatrix() {}

    double operator()(uint i, uint j) const THROW {
	check_indices(i,j);
	return data[i*ncol+j];
    }
    double& operator()(uint i, uint j) THROW {
	check_indices(i,j);
	return data[i*ncol+j];
    }

    bool operator==(const DMatrix& m) const THROW;
    bool operator!=(const DMatrix& m) const THROW;

    const DMatrix& operator+() const;
    DMatrix operator-() const;
    DMatrix operator+(const DMatrix& m) const THROW;
    DMatrix operator-(const DMatrix& m) const THROW;
    DMatrix operator*(double f) const;
    DMatrix operator/(double f) const THROW;
    DMatrix operator*(const DMatrix& m) const THROW;
    friend DMatrix operator*(double f, const DMatrix& p);

    // const DMatrix& operator=(const DMatrix& v);
    const DMatrix& operator+=(const DMatrix& v) THROW;
    const DMatrix& operator-=(const DMatrix& v) THROW;
    const DMatrix& operator*=(double f);
    const DMatrix& operator/=(double f) THROW;

    DMatrix t() const;
    bool is_symmetric() const;

    double norm_F() const;

    const double* as_vector() const {
	return &data[0];
    }
    double* as_vector() {
	return &data[0];
    }

    bool is_factored() const {
	return factored;
    }

    // DMatrix inv() const THROW;
    friend void dpotrf(DMatrix& A);
    friend void dposv(DMatrix& A, DMatrix& B);
    friend void dpotri(DMatrix& A);

    friend std::ostream& operator<<(std::ostream& os, const DMatrix& m);

    Vector operator*(const Vector& v) const THROW;
};

void dpotrf(DMatrix& A);
void dposv(DMatrix& A, DMatrix& B);

// void dgemm(double alpha, const DMatrix& A, CBLAS_TRANSPOSE opA, const DMatrix& B, CBLAS_TRANSPOSE opB,
	   // double beta, DMatrix& C);

class MapMatrix : public Matrix {
private:
    typedef std::map<uint,double> Row;
    std::vector<Row> data;

public:
    MapMatrix(uint n) {
	nrow = ncol = n;
	data.resize(n);
    }
    MapMatrix(const CSRMatrix& sm);

    double operator()(uint i, uint j) const THROW {
	check_indices(i,j);

	Row::const_iterator it = data[i].find(j);
	if (it == data[i].end())
	    return 0.0;
	return it->second;
    }
    double& operator()(uint i, uint j) THROW {
	check_indices(i,j);
	return data[i][j];
    }

    const std::map<uint,double>& operator()(uint i) const THROW {
	ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
	return data[i];
    }
    std::map<uint,double>& operator()(uint i) THROW {
	ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
	return data[i];
    }

    friend class CSRMatrix;
    friend class SkylineMatrix;

    friend std::ostream& operator<<(std::ostream& os, const MapMatrix& m);
};

class SkylineMatrix;
// We always assume that each row of CSRMatrix contains at least one element
class CSRMatrix : public Matrix {
protected:
    uvector<uint> ia, ja;
    uvector<double> a;

public:
    CSRMatrix();
    CSRMatrix(const MapMatrix& A);
    virtual ~CSRMatrix() { }

    const CSRMatrix& operator=(const CSRMatrix& A);
    uint size() const THROW {
	ASSERT(nrow == ncol, "size was called for rectangular matrix: nrow = " << nrow << ", ncol = " << ncol);
	return nrow;
    }

    virtual uint index(uint i, uint j) const;

    uvector<uint>& get_ia() { return ia; }
    uvector<uint>& get_ja() { return ja; }
    uvector<double>& get_a() { return a; }

    const uvector<uint>& get_ia() const { return ia; }
    const uvector<uint>& get_ja() const { return ja; }
    const uvector<double>& get_a() const { return a; }

    Vector operator*(const Vector& v) const THROW {
	Vector x(nrow);
	multiply(*this, v, x);
	return x;
    }

    void set_size(uint r, uint c) { nrow = r; ncol = c; }

    virtual double  operator()(uint i, uint j) const THROW;
    virtual double& operator()(uint i, uint j) THROW;
    virtual bool    exist(uint i, uint j) const THROW;

    virtual void load(const std::string& filename, DumpType type = BINARY) THROW;

    virtual bool is_symmetric() const;
    virtual std::string stat(bool ignore_pos_offdiagonal) const;

    double row_sum(uint i) const THROW;

    virtual uint nnz() const { return ia.back(); }

    virtual void reserve(uint n, uint nnz) {
	ia.reserve(n+1);
	ja.reserve(nnz);
	a.reserve(nnz);
    }

    // Friend classes
    friend class Prec;
    friend class SymPrec;
    friend class RelPrec;
    friend class AMGPrec;
    friend class SPEMesh;
    friend class MultiSplitPrec;

    // Friend functions
    friend void	multiply(const CSRMatrix& A, const Vector& v, Vector& res, char type = 'o') THROW;
    friend void transpose(const CSRMatrix& A, CSRMatrix& B);
    friend std::ostream& operator<<(std::ostream& os, const CSRMatrix& sm);

    friend void get_matrix_part(const SkylineMatrix& A, SkylineMatrix& lA, uint start, uint end);
    friend void construct_local_matrix(const SkylineMatrix& A, SkylineMatrix& locA,
				       std::vector<uint>& r, uint my_rank, uint ncpus);
};
void dump(const std::string& filename, const CSRMatrix& A, DumpType type) THROW;
void scale_c(CSRMatrix& A, double alpha);

/* Convert CSR matrix to CSC */
void convert(const CSRMatrix& A, int* ia, int* ja, double* a);

class SkylineMatrix : public CSRMatrix {
    uvector<uint> iasym, jasym;
    uvector<double> asym;

public:
    SkylineMatrix();
    SkylineMatrix(const MapMatrix& A);

    void optimize_storage(char type = 's');

    void load(const std::string& filename, DumpType type, bool transform) THROW;

    void permute(const std::vector<uint>& perm) THROW;

    uint index(uint i, uint j) const;

    friend void sym_multiply(const SkylineMatrix& A, const Vector& v, Vector& res) THROW;
};
void sym_multiply(const SkylineMatrix& A, const Vector& v, Vector& res) THROW;

void dump(const std::string& filename, const CSRMatrix& A, DumpType type) THROW;

/* r = b - Ax */
void residual(const CSRMatrix& A, const Vector& b, const Vector& x, Vector& r) THROW;

#endif // #ifndef __MATRIX_H__
