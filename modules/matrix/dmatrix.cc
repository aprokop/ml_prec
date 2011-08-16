#include "matrix.h"
#include "include/logger.h"
#include "include/tools.h"

#include <cmath>
#include <iostream>
#include <cstring>

using namespace std;

DEFINE_LOGGER("DMatrix");

DMatrix::DMatrix(uint n_, Type type) THROW {
    ASSERT(n_, "Can't create matrix 0x0");
    nrow = ncol = n_;
    data.resize(nrow*ncol);
    factored = false;

    switch (type) {
	case NONE: break;
	case IDENTITY: /* identity matrix */
	    memset(&data[0], 0, nrow*ncol*sizeof(double));
	    for (uint i = 0; i < nrow; i++)
		(*this)(i,i) = 1.0;
	    break;
	case ZERO:
	    memset(&data[0], 0, nrow*ncol*sizeof(double));
    }
}

DMatrix::DMatrix(uint m_, uint n_, Type type) THROW {
    ASSERT(m_ && n_, "Can't create matrix " << m_ << "x" << n_);
    nrow = m_;
    ncol = n_;
    data.resize(nrow*ncol);
    factored = false;

    switch (type) {
	case NONE: break;
	case IDENTITY: /* identity matrix */
	    ASSERT(m_ == n_, "IDENTITY is only for square matrices");
	    for (uint i = 0; i < nrow; i++)
		(*this)(i,i) = 1.0;
	    break;
	case ZERO:
	    memset(&data[0], 0, nrow*ncol*sizeof(double));
    }
}


DMatrix::DMatrix(const DMatrix& v) : Matrix() {
    nrow = v.nrow;
    ncol = v.ncol;
    data.resize(nrow*ncol);
    memcpy(&data[0], &v.data[0], nrow*ncol*sizeof(double));
    factored = false;
}

DMatrix::DMatrix(const char* values) THROW {
    ASSERT(values != NULL, "NULL string for constructor");

    std::istringstream buffer(values);
    double b;

    nrow = ncol = 0;
    while (buffer.peek() != EOF) {
	uint cols = 0;
	nrow++;
	while (buffer.peek() != EOF && buffer.peek() != ';') {
	    cols++;
	    buffer >> b;
	    data.push_back(b);
	    while (buffer.peek() == ' ')
		buffer.get();
	}
	if (ncol == 0)
	    ncol = cols;
	else {
	    ASSERT(ncol == cols, "Different number of elements in rows");
	}

	if (buffer.peek() == ';') {
	    buffer.get();
	    ASSERT(buffer.peek() != EOF, ";<EOF> in the input, wrong format");
	}
    }
    factored = false;
}

DMatrix::DMatrix(const std::vector<double>& v, uint _ncol) THROW {
    ASSERT(_ncol, "Trying to create matrix with 0 columns");
    ASSERT(v.size() % _ncol == 0, "Can not get integer number of rows: ncol = " << ncol << ", vector = " << v.size());
    ncol = _ncol;
    nrow = v.size() / ncol;
    data = v;
    factored = false;
}

DMatrix::DMatrix(const CSRMatrix& sm) {
    data.resize(sm.rows()*sm.cols());
    nrow = sm.rows();
    ncol = sm.cols();

    /* FIXME: this must be done much faster by going
     * through nonzero elements of CSRMatrix */
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    data[i*ncol+j] = sm(i,j);
    factored = false;
}

bool DMatrix::operator==(const DMatrix& m) const THROW {
    ASSERT(nrow == m.nrow && ncol == m.ncol, "Comparing matrices with different dimensions");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    return data == m.data;
}

bool DMatrix::operator!=(const DMatrix& m) const THROW {
    return !(*this == m);
}

const DMatrix& DMatrix::operator+() const {
    return *this;
}

DMatrix DMatrix::operator-() const {
    ASSERT(factored == false, "This operation is not available for factored matrices");
    DMatrix m(nrow, ncol);
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    m.data[i*ncol+j] = -data[i*ncol+j];
    return m;
}

DMatrix DMatrix::operator+(const DMatrix& m) const THROW {
    ASSERT(nrow == m.nrow && ncol == m.ncol, "Different dimensions of adding matrices");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    DMatrix r(nrow, ncol);
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    r.data[i*ncol+j] = data[i*ncol+j] + m.data[i*ncol+j];
    return r;
}

DMatrix DMatrix::operator-(const DMatrix& m) const THROW {
    ASSERT(nrow == m.nrow && ncol == m.ncol, "Different dimensions of adding matrices");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    DMatrix r(nrow, ncol);
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    r.data[i*ncol+j] = data[i*ncol+j] - m.data[i*ncol+j];
    return r;
}

DMatrix DMatrix::operator*(double f) const {
    ASSERT(factored == false, "This operation is not available for factored matrices");
    DMatrix r(nrow, ncol);
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    r.data[i*ncol+j] = f*data[i*ncol+j];
    return r;
}

DMatrix DMatrix::operator/(double f) const THROW {
    ASSERT(f != 0, "Division by zero");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    DMatrix r(nrow, ncol);
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    r.data[i*ncol+j] = data[i*ncol+j]/f;
    return r;
}

DMatrix DMatrix::operator*(const DMatrix& w) const THROW {
    ASSERT(ncol == w.nrow, "Not matched sizes: " << nrow << "x" << ncol << " and " << w.nrow << "x" << w.ncol);
    ASSERT(factored == false, "This operation is not available for factored matrices");
    int m = nrow, n = w.ncol, l = ncol;
    int i, j, k;
    double s;

    DMatrix r(m, n);
    for (i = 0; i < m; i++)
	for (j = 0; j < n; j++) {
	    for (k = 0, s = 0; k < l; k++) {
		s += data[i*l + k] * w.data[k*n + j];
	    }
	    r.data[i*n+j] = s;
	}

    return r;
}

// const DMatrix& DMatrix::operator=(const DMatrix& m) {
    // nrow = m.nrow;
    // ncol = m.ncol;
    // data = m.data;
    // return *this;
// }

const DMatrix& DMatrix::operator+=(const DMatrix& m) THROW {
    ASSERT(nrow == m.nrow && ncol == m.ncol, "Trying to add matrix with different dimensions");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    for(uint i = 0; i < nrow*ncol; i++)
	data[i] += m.data[i];
    return *this;
}

const DMatrix& DMatrix::operator-=(const DMatrix& m) THROW {
    ASSERT(nrow == m.nrow && ncol == m.ncol, "Trying to subtract matrix with different dimensions");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    for(uint i = 0; i < nrow*ncol; i++)
	data[i] -= m.data[i];
    return *this;
}

const DMatrix& DMatrix::operator*=(double f) {
    ASSERT(factored == false, "This operation is not available for factored matrices");
    for (uint i = 0; i < nrow*ncol; i++)
	data[i] *= f;
    return *this;
}

DMatrix operator*(double f, const DMatrix& p) {
    return p*f;
}

const DMatrix& DMatrix::operator/=(double f) THROW {
    ASSERT(is_not_equal(f,0), "Division by zero");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    for (uint i = 0; i < nrow*ncol; i++)
	data[i] /= f;
    return *this;
}

DMatrix DMatrix::t() const {
    ASSERT(factored == false, "This operation is not available for factored matrices");
    DMatrix r(ncol, nrow);
    for (uint i = 0; i < ncol; i++)
	for (uint j = 0; j < nrow; j++)
	    r.data[i*nrow+j] = data[j*ncol+i];
    return r;
}

std::ostream& operator<<(std::ostream& os, const DMatrix& w) {
    int m = w.rows(), n = w.cols();

    os << "Size: " << m << " x " << n << std::endl;
#if 0
    os << std::scientific;
#else
    os.precision(5);
    os.setf(std::ios::showpos);
    os << std::fixed;
#endif
    int i, j;
    if (w.is_factored() == false) {
	for (i = 0; i < m; i++, os << std::endl) {
	    os << "   ";
	    for (j = 0; j < n; j++)
		if (w(i,j))
		    os << w(i,j) << " ";
		else
		    os << "         ";
	}
    } else {
	for (i = 0; i < m; i++, os << std::endl) {
	    os << "   ";
	    for (j = 0; j < i; j++)
		os << w(j,i) << " ";
	    for (j = i; j < n; j++)
		os << w(i,j) << " ";
	}
	os << "Matrix is factored with the following factor:" << std::endl;
	for (i = 0; i < m; i++, os << std::endl) {
	    os << "   ";
	    for (j = 0; j < i; j++)
		os << w(i,j) << " ";
	    os << w.d[i];
	}
    }

    return os;
}

Vector DMatrix::operator*(const Vector& v) const THROW {
    ASSERT(ncol == v.size(), "Number of columns in matrix doesn't match vector dimension");
    ASSERT(factored == false, "This operation is not available for factored matrices");
    Vector r(nrow);
    double s;
    for (uint i = 0; i < nrow; i++) {
	s = 0.;
	for (uint j = 0; j < ncol; j++)
	    s += data[i*ncol+j] * v[j];
	r[i] = s;
    }
    return r;
}

double DMatrix::norm_F() const {
    ASSERT(factored == false, "This operation is not available for factored matrices");
    // TODO: use dnrm2 for this
    double s = 0.;
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < ncol; j++)
	    s += data[i*ncol+j]*data[i*ncol+j];
    return sqrt(s);
}

bool DMatrix::is_symmetric() const {
    ASSERT(nrow == ncol, "A is not square: " << sizes());
    for (uint i = 0; i < nrow; i++)
	for (uint j = 0; j < i; j++)
	    if (is_not_equal((*this)(i,j),(*this)(j,i)))
		return false;
    return true;
}

/*
 * Given an s.p.d matrix A routine constructs its Cholesky decomposition
 * A = L L^T. Only the upper triangle of A is used for input; it is not
 * modified. The Cholesky factor L is stored in the lower triangle of A,
 * except for its diagonal which is stored in d
 */
void dpotrf(DMatrix& A) {
    ASSERT(A.rows() == A.cols(), "A must be square: " << A.sizes());
    if (A.is_factored())
	return;
    uint n = A.rows();

    if (n == 1) {
	A.factored = true;
	return;
    }

    A.d.resize(n);
    for (uint i = 0; i < n; i++)
	for (uint j = i; j < n; j++) {
	    double sum = A(i,j);

	    for (int k = i-1; k >= 0; k--)
		sum -= A(i,k)*A(j,k);

	    if (i == j) {
		if (sum <= 0.0)
		    THROW_EXCEPTION("Failed: sum = " << sum);
		A.d[i] = sqrt(sum);
	    } else {
		A(j,i) = sum / A.d[i];
	    }
	}

    A.factored = true;
}

/*
 * Computes the solution of A * X = B with A being s.p.d
 * TODO: rewrite it w.r.t. daxpy and so on with proper increments
 */
void dposv(DMatrix& A, DMatrix& B) {
    ASSERT(A.rows() == A.cols(), "A must be square: " << A.sizes());
    ASSERT(B.rows() == A.rows(), "A: " << A.sizes() << ", B: " << B.sizes());
    uint n = A.rows(), nrhs = B.cols();
    if (A.is_factored() == false)
	dpotrf(A);

    if (n == 1) {
	double inv = 1./A(0,0);
	for (uint k = 0; k < nrhs; k++)
	    B(0,k) *= inv;
	return;
    }

    for (uint k = 0; k < nrhs; k++) {
	/* Solve system with L */
	for (uint i = 0; i < n; i++) {
	    for (uint j = 0; j < i; j++)
		B(i,k) -= A(i,j)*B(j,k);
	    B(i,k) /= A.d[i];
	}
	/* Solve system with L^T */
	for (int i = n-1; i >= 0; i--) {
	    for (int j = n-1; j > i; j--)
		/* Note the transposed indices in A */
		B(i,k) -= A(j,i)*B(j,k);
	    B(i,k) /= A.d[i];
	}
    }
}

void dpotri(DMatrix& A) {
    ASSERT(A.rows() == A.cols(), "A must be square: " << A.sizes());

    uint n = A.rows();
    DMatrix B(n, DMatrix::IDENTITY), TMP(A);

    dposv(A,B); /* B <- A^{-1} */
    memcpy(&A.data[0], &B.data[0], n*n*sizeof(double));
    A.factored = false;
}

// #ifdef HAVE_LAPACK
// void dposv(DMatrix& A, DMatrix& B) {
    // ASSERT(A.rows() == A.cols(), A.sizes());
    // ASSERT(B.rows() == A.rows(), "A: " << A.sizes() << ", B: " << B.sizes());
//
    // std::cout << A << std::endl;
    // std::cout << B << std::endl;
//
    // // int clapack_dposv(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
		      // // const int N, const int NRHS, double *A, const int lda,
		      // // double *B, const int ldb);
    // // int info = clapack_dposv(CblasRowMajor, CblasUpper, m, n, A.as_vector(), A.rows(), B.as_vector(), B.rows());
    // // int info = clapack_dposv(CblasRowMajor, CblasUpper, A.rows(), B.cols(), A.as_vector(), A.rows(), B.as_vector(), B.rows());
    // // int info = clapack_dposv(CblasColMajor, CblasUpper, A.rows(), B.cols(), A.as_vector(), A.rows(), B.as_vector(), B.rows());
    // int info = clapack_dposv(CblasColMajor, CblasUpper, A.rows(), B.cols(), A.as_vector(), A.rows(), B.as_vector(), B.cols());
    // ASSERT(info == 0, "lapack dposv returned " << info);
//
    // std::cout << B << std::endl;
    // exit(1);
// }
//
// void dgemm(double alpha, const DMatrix& A, CBLAS_TRANSPOSE opA, const DMatrix& B, CBLAS_TRANSPOSE opB,
	   // double beta, DMatrix& C) {
    // uint m = (opA == CblasNoTrans) ? A.rows() : A.cols();
    // uint n = (opB == CblasNoTrans) ? B.cols() : B.rows();
    // uint k = (opA == CblasNoTrans) ? A.cols() : A.rows();
    // ASSERT(C.rows() == m && C.cols() == n, "C: " << C.sizes() << ", op(A)op(B): " << m << "x" << n);
//
    // cblas_dgemm(CblasRowMajor, opA, opB, m, n, k, alpha,
		// A.as_vector(), A.rows(), B.as_vector(), B.rows(), beta, C.as_vector(), C.rows());
// }
// #else
// #error "NOT_IMPLEMENTED"
// #endif
