#include "matrix.h"
#include "include/logger.h"
#include "include/tools.h"

DEFINE_LOGGER("MapMatrix");

MapMatrix::MapMatrix(const CSRMatrix& sm) {
    const uint n = sm.size();
    const uvector<uint> &ia = sm.get_ia(), &ja = sm.get_ja();
    const uvector<double> &a = sm.get_a();

    data.resize(n);
    for (uint i = 0; i < n; i++)
	for (uint j = ia[i]; j < ia[i+1]; j++)
	    data[i][ja[j]] = a[j];
}

void multiply(const MapMatrix& A, const Vector& v, Vector& res) THROW {
    uint n = A.rows();
    ASSERT(res.size() == A.nrow, "Not enough space in res vector");
    if (res.size() == 0)
	return;

    ASSERT(A.rows() == res.size(), "Different sizes: A is " << A.sizes() << ", res is " << res.size());
    ASSERT(A.ncol == v.size(), "Multiplying sparse matrix and vector with different dimensions");

    memset(&res[0], 0, res.size()*sizeof(double));
    for (uint i = 0; i < A.nrow; i++) {
	const MapMatrix::Row& row = A.data[i];
	for (MapMatrix::Row::const_iterator it = row.begin(); it != row.end(); it++)
	    res[i] += it->second * v[it->first];
    }
}

std::ostream& operator<<(std::ostream& os, const MapMatrix& sm) {
    os << "Size: " << sm.nrow << "x" << sm.ncol << std::endl;
    for (uint i = 0; i < sm.nrow; i++) {
	os << "  Row: " << i << std::endl;
	const MapMatrix::Row& row = sm.data[i];
	for (MapMatrix::Row::const_iterator it = row.begin(); it != row.end(); it++)
	    os << "    " << it->first << ": " << it->second << std::endl;
	os << std::endl;
    }

    return os;
}
