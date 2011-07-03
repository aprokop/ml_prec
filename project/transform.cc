#include "modules/matrix/matrix.h"
#include "include/tools.h"

void transform(CSRMatrix& A_, Vector& b) {
    uint n = A_.size();
    MapMatrix A(A_), B(n);

    double alpha = 1.0;

    typedef std::map<uint,double> map;

    /* Scale matrix A so that A(k,k) = 1 */
    for (uint i = 0; i < n; i++) {
	double d = 1./A(i,i);
	for (map::iterator it = A(i).begin(); it != A(i).end(); it++)
	    it->second *= d;
    }

    for (uint i = 0; i < n; i++) {
	for (map::const_iterator it = A(i).begin(); it != A(i).end(); it++) {
	    uint k = it->first;
	    for (map::const_iterator it1 = A(k).begin(); it1 != A(k).end(); it1++)
		if (it->first != i) B(i,it1->first) += alpha*(-it->second)*it1->second;
		else		    B(i,it1->first) += it1->second;
	}
    }

    A_ = SkylineMatrix(B);
}
