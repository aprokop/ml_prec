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

    void check_indices(uint i, uint j) const{
        ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
        ASSERT(j < ncol, "Col index is out of boudaries: j = " << j << ", ncol = " << ncol);
    }

public:
    virtual ~Matrix() {}
    // virtual Vector operator*(const Vector& v) const = 0;

    uint rows() const { return nrow; }
    uint cols() const { return ncol; }
    uint size() const{
        ASSERT(nrow == ncol, "size was called for rectangular matrix: nrow = " << nrow << ", ncol = " << ncol);
        return nrow;
    }

    virtual double  operator()(uint i, uint j) const = 0;
    virtual double& operator()(uint i, uint j) = 0;

    std::pair<uint,uint> sizes() const {
        return std::pair<uint,uint>(nrow,ncol);
    }
};

class CSRMatrix;
/* Dense matrix class */
class DMatrix : public Matrix {
private:
    uvector<double> data;

    /*
     * Symmetric matrix can also be stored in factored form.
     * In this case the original matrix is stored in upper triangle,
     * L (except the diagonal) is stored in lower triangle, diagonal
     * of L is stored in d
     */
    bool factored;
    uvector<double> d;

public:
    enum Type {
        NONE,
        ZERO,
        IDENTITY
    };

public:
    DMatrix(uint n_ = 1, Type type = ZERO);
    DMatrix(uint m_, uint n_, Type type = NONE);
    DMatrix(const DMatrix& v);
    DMatrix(const char* values);
    DMatrix(const uvector<double>& v, uint _ncol);
    DMatrix(const CSRMatrix& sm);
    ~DMatrix() {}

    double operator()(uint i, uint j) const{
        check_indices(i,j);
        return data[i*ncol+j];
    }
    double& operator()(uint i, uint j){
        check_indices(i,j);
        return data[i*ncol+j];
    }

    bool operator==(const DMatrix& m) const;
    bool operator!=(const DMatrix& m) const;

    const DMatrix& operator+() const;
    DMatrix operator-() const;
    DMatrix operator+(const DMatrix& m) const;
    DMatrix operator-(const DMatrix& m) const;
    DMatrix operator*(double f) const;
    DMatrix operator/(double f) const;
    DMatrix operator*(const DMatrix& m) const;
    friend DMatrix operator*(double f, const DMatrix& p);

    // const DMatrix& operator=(const DMatrix& v);
    const DMatrix& operator+=(const DMatrix& v);
    const DMatrix& operator-=(const DMatrix& v);
    const DMatrix& operator*=(double f);
    const DMatrix& operator/=(double f);

    void get_submatrix(const uvector<uint>& r, const uvector<uint>& c, DMatrix& D1) const;

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

    DMatrix inv() const;
    friend void dpotrf(DMatrix& A);
    friend void dposv(DMatrix& A, DMatrix& B);
    friend void dpotri(DMatrix& A);

    friend void	multiply(const DMatrix& A, const Vector& v, Vector& res);
    friend std::ostream& operator<<(std::ostream& os, const DMatrix& m);

    Vector operator*(const Vector& v) const;
};

void dpotrf(DMatrix& A);
void dposv(DMatrix& A, DMatrix& B);

/* Compute eigenvalues of a s.p.d 3x3 matrix */
void eigs3(const DMatrix& A, double r[3]);


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

    void set_size(uint r, uint c) {
        nrow = r;
        ncol = c;
        data.resize(r);
    }

    double operator()(uint i, uint j) const{
        check_indices(i,j);

        Row::const_iterator it = data[i].find(j);
        if (it == data[i].end())
            return 0.0;
        return it->second;
    }
    double& operator()(uint i, uint j){
        check_indices(i,j);
        return data[i][j];
    }

    const std::map<uint,double>& operator()(uint i) const{
        ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
        return data[i];
    }
    std::map<uint,double>& operator()(uint i){
        ASSERT(i < nrow, "Row index is out of boudaries: i = " << i << ", nrow = " << nrow);
        return data[i];
    }

    friend class CSRMatrix;
    friend class SkylineMatrix;

    friend void	multiply(const MapMatrix& A, const Vector& v, Vector& res);
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
    CSRMatrix(const DMatrix& A);
    virtual ~CSRMatrix() { }

    void set_size(uint r, uint c) { nrow = r; ncol = c; }

    const CSRMatrix& operator=(const CSRMatrix& A);

    virtual uint index(uint i, uint j) const;

    uvector<uint>& get_ia() { return ia; }
    uvector<uint>& get_ja() { return ja; }
    uvector<double>& get_a() { return a; }

    const uvector<uint>& get_ia() const { return ia; }
    const uvector<uint>& get_ja() const { return ja; }
    const uvector<double>& get_a() const { return a; }

    Vector operator*(const Vector& v) const{
        Vector x(nrow);
        multiply(*this, v, x);
        return x;
    }

    virtual double  operator()(uint i, uint j) const;
    virtual double& operator()(uint i, uint j);
    virtual bool    exist(uint i, uint j) const;

    virtual void load(const std::string& filename, DumpType type = BINARY);

    virtual bool is_symmetric() const;
    virtual std::string stat(bool ignore_pos_offdiagonal) const;

    double row_sum(uint i) const;

    virtual uint nnz() const { return ((nrow && ncol) ? ia.back() : 0); }

    void swap(CSRMatrix& B) {
        ia.swap(B.ia);
        ja.swap(B.ja);
        a.swap(B.a);
        std::swap(nrow, B.nrow);
        std::swap(ncol, B.ncol);
    }

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
    friend class Mesh;
    friend class MultiSplitPrec;

    // Friend functions
    friend void	multiply(const CSRMatrix& A, const Vector& v, Vector& res, char type = 'o');
    friend void transpose(const CSRMatrix& A, CSRMatrix& B);
    friend std::ostream& operator<<(std::ostream& os, const CSRMatrix& sm);

    friend void get_matrix_part(const SkylineMatrix& A, SkylineMatrix& lA, uint start, uint end);
    friend void construct_local_matrix(const SkylineMatrix& A, SkylineMatrix& locA,
                                       std::vector<uint>& r, uint my_rank, uint ncpus);
};
void dump(const std::string& filename, const CSRMatrix& A, DumpType type);
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

    void load(const std::string& filename, DumpType type, bool transform);

    void permute(const std::vector<uint>& perm);

    uint index(uint i, uint j) const;

    void swap(SkylineMatrix& B) {
        CSRMatrix::swap(B);
        iasym.swap(B.iasym);
        jasym.swap(B.jasym);
    }
    friend void sym_multiply(const SkylineMatrix& A, const Vector& v, Vector& res);
};
void sym_multiply(const SkylineMatrix& A, const Vector& v, Vector& res);

void dump(const std::string& filename, const CSRMatrix& A, DumpType type);

/* r = b - Ax */
void residual(const CSRMatrix& A, const Vector& b, const Vector& x, Vector& r);

#endif // #ifndef __MATRIX_H__
