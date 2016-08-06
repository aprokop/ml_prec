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

DMatrix::DMatrix(const uvector<double>& v, uint ncol_) THROW {
    ASSERT(ncol_, "Trying to create matrix with 0 columns");
    ASSERT(v.size() % ncol_ == 0, "Can not get integer number of rows: ncol = " << ncol_ << ", vector = " << v.size());
    ncol = ncol_;
    nrow = v.size() / ncol;
    data = v;
    factored = false;
}

DMatrix::DMatrix(const CSRMatrix& sm) {
    data.resize(sm.rows()*sm.cols());
    nrow = sm.rows();
    ncol = sm.cols();

    const uvector<uint>& ia = sm.get_ia();
    const uvector<uint>& ja = sm.get_ja();
    const uvector<double>& a = sm.get_a();
    for (uint i = 0; i < nrow; i++)
	for (uint j = ia[i]; j < ia[i+1]; j++)
	    data[i*ncol + ja[j]] = a[j];
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
    uint m = nrow, n = w.ncol, l = ncol;
    uint i, j, k;
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

/* Arrays r and c may be unsorted, but each index must be unique */
void DMatrix::get_submatrix(const uvector<uint>& r, const uvector<uint>& c, DMatrix& D1) const {
    D1 = DMatrix(r.size(), c.size());
    for (uint i = 0; i < r.size(); i++)
	for (uint j = 0; j < c.size(); j++)
	    D1(i,j) = (*this)(r[i], c[j]);
}

/* Horrowfull implementation. For n > 10 is not useful */
DMatrix DMatrix::inv() const THROW {
    ASSERT(nrow == ncol, "Matrix is not square: " << nrow << " x " << ncol);
    uint n = nrow;

    if (n == 1) {
	if (is_equal(data[0], 0.0))
	    throw Exception("Inverting singular matrix");

	DMatrix I(1);
	I(0,0) = 1/data[0];
	return I;
    }

    DMatrix R(*this);
    DMatrix Q(n, DMatrix::IDENTITY);
    DMatrix U(n);
    std::vector<double> x(n);

    /* Find A = QR using reflection method */
    for (uint i = 0; i < n-1; i++) {
	x[i] = R(i,i);

	double s = 0.0;
	for (uint j = i+1; j < n; j++) {
	    x[j] = R(j,i);
	    s += x[j]*x[j];
	}

	double norm = sqrt(s + x[i]*x[i]);
	if (fabs(norm < 1e-14))
	    throw Exception("Inverting singular matrix");
	if (is_equal(s, 0.0))
	    continue;

	x[i] -= norm;
	norm = sqrt(s + x[i]*x[i]);

	for (uint j = i; j < n; j++)
	    x[j] /= norm;

	/* Generate U matrix */
	for (uint k = 0; k < n; k++)
	    for (uint l = 0; l < n; l++)
		if (k < i || l < i)
		    U(k,l) = (k == l) ? 1 : 0;
		else
		    U(k,l) = (k == l) ? 1 - 2*x[k]*x[l] : -2*x[k]*x[l];

	/* Modify Q matrix */
	Q = Q * U;
	/* Modify R matrix */
	R = U*R;
    }
    /* We need to check last element for non singularity */
    if (fabs(R(n-1,n-1)) < 1e-14)
	throw Exception("Inversing singular matrix");

#if 0
    {
	/* Check */
	double tmp, max = 0;
	DMatrix QR = Q*R;
	for (uint i = 0; i < n; i++)
	    for (uint j = 0; j < n; j++)
		if ((tmp = fabs(QR.data[i*n+j] - data[i*n+j])) > max)
		    max = tmp;
	LOG_INFO("max(QR - A) : " << max);
    }
#endif

    /* Find the inverse R matrix. It is done simply through (R | E) algorithm */
    DMatrix R1(n, DMatrix::IDENTITY);
    for (uint i = 0; i < n; i++)
	R1(i,i) /= R(i,i);
    for (int i = n-2; i >= 0; i--)
	/* Subtract i+1,i+2,...,n-1 lines from i-th */
	for (uint k = i+1; k < n; k++)
	    for (uint j = k; j < n; j++)
		R1(i,j) -= (R(i,k)/R(i,i))*R1(k,j);

#if 0
    {
	/* Check the inverse R */
	double tmp, max = 0;
	DMatrix I = R*R1;
	for (uint i = 0; i < n; i++)
	    for (uint j = 0; j < n; j++)
		if ((tmp = fabs(I(i,j) - ((i == j) ? 1 : 0))) > max)
		    max = tmp;
	LOG_INFO("max(R*R1 - I) : " << max);
    }
#endif

    return R1*Q.t();
}

std::ostream& operator<<(std::ostream& os, const DMatrix& w) {
    uint m = w.rows(), n = w.cols();

    os << "Size: " << m << " x " << n << std::endl;
#if 1
    os << std::scientific;
#else
    os.precision(5);
    os.setf(std::ios::showpos);
    os << std::fixed;
#endif
    uint i, j;
    if (w.is_factored() == false) {
	for (i = 0; i < m; i++, os << std::endl) {
	    os << "   ";
	    for (j = 0; j < n; j++)
#if 0
		if (w(i,j))
		    os << w(i,j) << " ";
		else
		    os << "         ";
#else
		os << w(i,j) << " ";
#endif
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
    os.unsetf(std::ios::showpos);

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

void multiply(const DMatrix& A, const Vector& v, Vector& res) THROW {
    ASSERT(A.rows() == res.size(), "Different sizes: A is " << A.sizes() << ", res is " << res.size());
    if (res.size() == 0)
	return;

    uint m = A.nrow, n = A.ncol;
    for (uint i = 0; i < m; i++) {
	res[i] = 0.0;
	for (uint j = 0; j < n; j++)
	    res[i] += A.data[i*n+j]*v[j];
    }
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
		sum -= A(i,uint(k))*A(j,uint(k));

	    if (i == j) {
		if (sum <= 0.0)
		    THROW_EXCEPTION("Failed: sum = " << sum << " at row " << i);
		A.d[i] = sqrt(sum);
	    } else {
		A(j,i) = sum / A.d[i];
	    }
	}

    A.factored = true;
}

/*
 * Computes the solution of A * X = B with A being s.p.d, storing X in place of B
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

/*
 * To find roots of the following cubic equation, where a, b, c and d are real
 *	a*x^3 + b*x^2 + c*x + d = 0
 */
static void cubic_roots(double a, double b, double c, double d, double r[3], bool all_are_real = true) {
    double p, q;
    double D;
    double u, v;

    /* Calculate p and q */
    p = (3*c/a - pow(b/a,2)) / 3;
    q = (2*pow(b/a,3) - 9*b*c/a/a + 27*d/a) / 27;

    /* Calculate discriminant D */
    D = pow(p/3,3) + pow(q/2,2);

    if (D > 0 && all_are_real == false) {
	/* D > 0 : one real and two complex roots */
	THROW_EXCEPTION("Complex roots: D = " << D);

	/* Calculate u and v */
	u = cbrt(-q/2 + sqrt(D));
	v = cbrt(-q/2 - sqrt(D));

	/* Find the three transformed roots */
	// r[0] = u + v;
	// r[1] = -(u + v)/2 + i (u-v)*sqrt(3)/2;
	// r[2] = -(u + v)/2 - i (u-v)*sqrt(3)/2;
    } else if (D > -1e-14) {
	/* D = 0 : three real roots, of which at least two are equal */

	/* Calculate u and v */
	u = v = cbrt(-q/2);

	/* Find the three transformed roots */
	r[0] = u + v;
	r[1] = r[2] = -(u + v)/2;
    } else {
	/* D < 0 : three distinct real roots */
	p = fabs(p)/3;
	double phi = acos(-q/2/sqrt(pow(p,3)));

        r[0] =  2 * sqrt(p) * cos(phi/3);
        r[1] = -2 * sqrt(p) * cos((phi+PI)/3);
        r[2] = -2 * sqrt(p) * cos((phi-PI)/3);
    }

    for (uint i = 0; i < 3; i++)
	r[i] -= b/a/3;
}

/* Compute eigenvalues of a s.p.d 3x3 matrix */
void eigs3(const DMatrix& A, double r[3]) {
    /*
     * Compute coefficients of a characteristic polynomial
     *	    \lambda^3 + A*\lambda^2 + B*\lambda + C = 0
     */
    double a, b, c;
    a = -(A(0,0) + A(1,1) + A(2,2));
    b = A(0,0)*A(1,1) + A(0,0)*A(2,2) + A(1,1)*A(2,2) - pow(A(0,1),2) - pow(A(0,2),2) - pow(A(1,2),2);
    c = -(A(0,0)*A(1,1)*A(2,2) + 2*A(0,1)*A(0,2)*A(1,2)) +
	    pow(A(0,1),2)*A(2,2) + pow(A(0,2),2)*A(1,1) + pow(A(1,2),2)*A(0,0);

    // LOG_VAR(N);
    // LOG_DEBUG(std::showpos << "x^3 " << A << "*x^2 " << B << "*x " << C << " = 0");

    cubic_roots(1., a, b, c, r);

    // LOG_DEBUG("Roots: " << r[0] << " " << r[1] << " " << r[2]);
}
