#include "matrix.h"
#include "include/logger.h"
#include "include/exception.h"
#include "include/tools.h"
#include "include/time.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <limits>
#include <typeinfo>

DEFINE_LOGGER("CSRMatrix");

CSRMatrix::CSRMatrix() {
    nrow = ncol = 0;
}

CSRMatrix::CSRMatrix(const MapMatrix& A) {
    nrow = A.rows();
    ncol = A.cols();
    ia.resize(nrow + 1);
    ja.reserve(nrow*7);
    a.reserve(nrow*7);

    ia[0] = 0;
    for (uint i = 0; i < nrow; i++) {
	const MapMatrix::Row& row = A.data[i];

	/* Add remaining elements of a row */
	for (MapMatrix::Row::const_iterator it = row.begin(); it != row.end(); it++) {
	    ja.push_back(it->first);
	    a.push_back(it->second);
	}

	ia[i+1] = ia[i] + row.size();
    }
}

CSRMatrix::CSRMatrix(const DMatrix& A) {
    nrow = A.rows();
    ncol = A.cols();
    ia.resize(nrow+1);
    ja.reserve(nrow*7);
    a.reserve(nrow*7);

    ia[0] = 0;
    for (uint i = 0; i < nrow; i++) {
	for (uint j = 0; j < ncol; j++)
	    if (A(i,j)) {
		ja.push_back(j);
		a.push_back(A(i,j));
	    }

	ia[i+1] = ja.size();
    }
}

const CSRMatrix& CSRMatrix::operator=(const CSRMatrix& A) {
    nrow = A.nrow;
    ncol = A.ncol;

    if (nrow) {
	uint n = A.ia.size();
	uint nnz = A.ia[n-1];
	ia.resize(n);
	ja.resize(nnz);
	a.resize(nnz);
	memcpy(ia.data(), A.ia.data(), n*sizeof(int));
	memcpy(ja.data(), A.ja.data(), nnz*sizeof(int));
	memcpy(a.data(),  A.a.data(),  nnz*sizeof(double));
    } else {
	ia.clear();
	ja.clear();
	a.clear();
    }

    return *this;
}

uint CSRMatrix::index(uint i, uint j) const {
    check_indices(i,j);

    uvector<uint>::const_iterator start = ja.begin() + ia[i];
    uvector<uint>::const_iterator   end = ja.begin() + ia[i+1];
    uvector<uint>::const_iterator it = std::lower_bound(start, end, j);
    if (it != end && !(j < *it))
	return it-ja.begin();

    return uint(-1);
}

bool CSRMatrix::exist(uint i, uint j) const THROW {
    return index(i,j) != uint(-1);
}

double CSRMatrix::operator()(uint i, uint j) const THROW {
    uint ind = index(i,j);
    if (ind != uint(-1))
	return a[ind];

    // LOG_WARN("(" << i << "," << j << ") is not in stencil");
    return 0;
}

double& CSRMatrix::operator()(uint i, uint j) THROW {
    uint ind = index(i,j);
    ASSERT(ind != uint(-1), "(" << i << "," << j << ") is not in stencil");

    return a[ind];
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

double CSRMatrix::row_sum(uint i) const THROW {
    ASSERT(nrow*ncol, "One of matrix dimensions is zero");
    return std::accumulate(&a[ia[i]], &a[ia[i+1]], 0.0);
}

bool CSRMatrix::is_symmetric() const {
    ASSERT(nrow == ncol, "Matrix is not square: " << sizes());
    for (uint i = 0; i < nrow; i++)
	for (uint j = ia[i]; j < ia[i+1]; j++)
	    if (a[j] != (*this)(ja[j], i))
		return false;
    return true;
}

void CSRMatrix::load(const std::string& filename, DumpType type) THROW {
    uint nnz;
    if (type == BINARY) {
	std::ifstream is(filename.c_str(), std::ifstream::binary);

	if (!is.good())
	    THROW_EXCEPTION("Problem reading file \"" << filename  << "\"");

	is.read(reinterpret_cast<char*>(&nrow), sizeof(uint));
	is.read(reinterpret_cast<char*>(&ncol), sizeof(uint));

	ia.resize(nrow+1);
	is.read(reinterpret_cast<char*>(&ia[0]), (nrow+1)*sizeof(uint));

	nnz = ia[nrow];
	ja.resize(nnz);
	a.resize(nnz);
	is.read(reinterpret_cast<char*>(&ja[0]), nnz*sizeof(uint));
	is.read(reinterpret_cast<char*>(&a[0]),  nnz*sizeof(double));
    } else if (type == ASCII) {
	std::ifstream is(filename.c_str());
	ASSERT(is.good(), "Problem reading file \"" << filename  << "\"");

	const int N = 2009;
	char str[N];

	is.getline(str, N); /* "# rows cols nonzeros" */
	is >> nrow >> ncol >> nnz;

	is.getline(str, N);
	is.getline(str, N); /* "# row_ptr" */
	ia.resize(nrow+1);
	for (uint i = 0; i <= nrow; i++)
	    is >> ia[i];

	is.getline(str, N);
	is.getline(str, N); /* "# col_ind" */
	ja.resize(nnz);
	for (uint i = 0; i < nnz; i++)
	    is >> ja[i];

	is.getline(str, N);
	is.getline(str, N); /* "# value" */
	a.resize(nnz);
	for (uint i  = 0; i < nnz; i++)
	    is >> a[i];
    } else if (type == MATRIX_MARKET) {
	std::ifstream is(filename.c_str());
	ASSERT(is.good(), "Problem reading file \"" << filename  << "\"");

	const int N = 2011;
	char str[N];
	is.getline(str, N); /* "%%MatrixMarket matrix coordinate real general" */
	/* TODO: ignore comments */
	is.getline(str, N); /* some comments */

	is >> nrow >> ncol >> nnz;

	MapMatrix mm(nrow);

	uint i, j;
	double v;
	for (uint k = 0; k < nnz; k++) {
	    is >> i >> j >> v;

	    mm(i-1,j-1) = v;
	}

	*this = CSRMatrix(mm);
    } else if (type == HYPRE) {
	THROW_EXCEPTION("Loading format HYPRE is not applicable to CSRMatrix");
    }

    LOG_INFO("Loaded matrix: sizes = " << sizes() << ", nnz = " << nnz);
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
    ASSERT(A.rows() == res.size(), "Different sizes: A is " << A.sizes() << ", res is " << res.size());
    if (res.size() == 0)
	return;

    memset(&res[0], 0, res.size()*sizeof(double));
    switch (type) {
	case 'o': /* A*v */
	    ASSERT(A.ncol == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.nrow, "Not enough space in res vector");

	    for (uint i = 0; i < A.nrow; i++)
		for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
		    res[i] += A.a[j] * v[A.ja[j]];
	    break;

	case 't': /* A^T*v */
	    ASSERT(A.nrow == v.size(), "Multiplying sparse matrix and vector with different dimensions");
	    ASSERT(res.size() == A.ncol, "Not enough space in res vector");

	    for (uint i = 0; i < A.nrow; i++)
		for (uint j = A.ia[i]; j < A.ia[i+1]; j++)
		    res[A.ja[j]] += A.a[j] * v[i];
	    break;
	case 's': /* A*v with A = A^T */
	    sym_multiply(dynamic_cast<const SkylineMatrix&>(A), v, res);
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

void dump(const std::string& filename, const CSRMatrix& A, DumpType type) THROW {
    const uvector<uint> &ia = A.get_ia(), &ja = A.get_ja();
    const uvector<double>& a = A.get_a();
    uint m = A.rows(), n = A.cols(), nnz = ia.back();

    std::ofstream os;
    if (type != BINARY) os.open(filename.c_str());
    else		os.open(filename.c_str(), std::ofstream::binary);

    if (!os.good())
	THROW_EXCEPTION("Could not open \"" << filename << "\" for write");

    if (type == BINARY) {
	os.write(reinterpret_cast<const char*>(&m), sizeof(uint));
	os.write(reinterpret_cast<const char*>(&n), sizeof(uint));
	os.write(reinterpret_cast<const char*>(&ia[0]), (m+1)*sizeof(uint));

	os.write(reinterpret_cast<const char*>(&ja[0]), nnz*sizeof(uint));
	os.write(reinterpret_cast<const char*>(&a[0]),  nnz*sizeof(double));
    } else if (type == ASCII) {
	os << "# rows cols nonzeros" << std::endl;
	os << m << " " << n << " " << nnz << std::endl;
	os << "# row_ptr" << std::endl;
	for (uint i = 0; i <= m; i++)
	    os << ia[i] << std::endl;
	os << "# col_ind" << std::endl;
	for (uint j = 0; j < nnz; j++)
	    os << ja[j] << std::endl;
	os << "# value" << std::endl;
	os << std::scientific << std::setprecision(15);
	for (uint j = 0; j < nnz; j++)
	    os << a[j] << std::endl;
    } else if (type == HYPRE) {
	try {
	    const SkylineMatrix& A_ = dynamic_cast<const SkylineMatrix&>(A);
	} catch (std::bad_cast) {
	    THROW_EXCEPTION("Dumping matrix in HYPRE format is not supported for CSRMatrix");
	}

	os << "0 " << m-1 << " 0 " << n-1 << std::endl;
	os << std::scientific << std::setprecision(15);
	for (uint i = 0; i < m; i++)
	    for (uint j = ia[i]; j < ia[i+1]; j++)
		os << i << " " << ja[j] << " " << a[j] << std::endl;
    } else if (type == MATRIX_MARKET) {
	os << std::fixed << std::setprecision(15);

	os << "%%MatrixMarket matrix coordinate real general" << std::endl;
	os << m << " " << n << " " << nnz << std::endl;
	for (uint i = 0; i < m; i++)
	    for (uint j = ia[i]; j < ia[i+1]; j++)
		os << i+1 << " " << ja[j]+1 << " " << a[j] << std::endl;
    }
}

void scale_c(CSRMatrix& A, double alpha) {
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    uvector<double>& a = A.get_a();
    uint n = A.rows();

    for (uint i = 0; i < n; i++) {
	double c = std::accumulate(&a[ia[i]], &a[ia[i+1]], 0.0);
	A(i,i) += (alpha-1)*c;
    }
}

std::string CSRMatrix::stat(bool ignore_pos_offdiagonal) const {
    std::ostringstream os;
    os << std::scientific;
#if 0
    const int N = 16;
    double ticks[N] = {-1e6,-1e3, -100, -10, -1, -0.1 -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1, 10, 100, 1e3, 1e6};
    int bins[N] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double s, pos = 0, neg = 0;
    for (uint i = 0; i < nrow; i++) {
	s = 0;
	for (uint j = ia[i]; j < ia[i+1]; j++) {
	    if (ja[j] != i) {
		// off-diagonal element
		bins[std::lower_bound(ticks, ticks+N, a[j]) - ticks]++;

		if (a[j] > 0) {
		    if (ignore_pos_offdiagonal == false)
			LOG_DEBUG("Positive element: (" << i << "," << ja[j] << ") : " << a[j]);
		    pos++;
		} else if (a[j] < 0) {
		    neg++;
		}
		s += a[j];
	    } else {
		if (a[j] < 0)
		    os << "Negative diagonal element: " << i << " : " << a[j] << std::endl;
		s += a[j];
	    }
	}
	if (s < -1e-8)
	    os << "Negative row sum: " << i << " : " << s << std::endl;
    }
    os << 100*pos/(pos+neg) << "% off-diagonal are positive" << std::endl;
    for (int i = 0; i < N-1; i++)
	os << "[" << ticks[i] << "," << ticks[i+1] << "] : " << bins[i] << std::endl;
#endif

    return os.str();
}

/* Convert CSR matrix to CSC */
void convert(const CSRMatrix& A, int* ia, int* ja, double* a) {
    const uint n = A.size();

    const uvector<uint>&   Aia = A.get_ia();
    const uvector<uint>&   Aja = A.get_ja();
    const uvector<double>&  Aa = A.get_a();

    std::vector<int> marker(n, 0);

    /* Calculate number of elements in each column */
    for (uint i = 0; i < n; i++)
	for (uint j = Aia[i]; j < Aia[i+1]; j++)
	    marker[Aja[j]]++;

    /*
     * Fill array ia for transposed matrix using marker
     * After cycle marker = B.ia (but shorter by 1 element)
     */
    ia[0] = 0;
    for (uint i = 0; i < n; i++) {
	ia[i+1] = ia[i] + marker[i];
	marker[i] = ia[i];
    }

    /* Fill in ja and a */
    int col, relpos;
    for (uint i = 0; i < n; i++)
	for (uint j = Aia[i]; j < Aia[i+1]; j++) {
	    col = Aja[j];
	    relpos = marker[col];

	    ja[relpos] = i;
	    a[relpos] = Aa[j];

	    marker[col]++;
	}
}
