#include "modules/matrix/matrix.h"
#include "include/tools.h"
#include "config.h"

void transform(CSRMatrix& A_, Vector& b, TransType transform) {
    uint n = A_.size();

    if (transform == TRANS_ILU || transform == TRANS_IL || transform == TRANS_IU) {
	MapMatrix A(A_), B(n);

	double alpha = 1.0;

	typedef std::map<uint,double> map;

	/* Scale matrix A so that A(k,k) = 1 */
	for (uint i = 0; i < n; i++) {
	    double d = 1./A(i,i);
	    for (map::iterator it = A(i).begin(); it != A(i).end(); it++)
		it->second *= d;
	}

	if (transform == TRANS_ILU) {
	    /* A -> (I + alpha*L + alpha*U)*A */
	    for (uint i = 0; i < n; i++) {
		for (map::const_iterator it = A(i).begin(); it != A(i).end(); it++) {
		    uint   k = it->first;
		    double v = (k != i) ? alpha*(-it->second) : 1.0;
		    for (map::const_iterator it1 = A(k).begin(); it1 != A(k).end(); it1++)
			B(i,it1->first) += v*it1->second;
		}
	    }
	}
	if (transform == TRANS_IL) {
	    /* A -> (I + alpha*L)*A */
	    alpha = 0.5;
	    for (uint i = 0; i < n; i++) {
		for (map::const_iterator it = A(i).begin(); it != A(i).end(); it++) {
		    uint k = it->first;
		    if (k <= i) {
			double v = (k != i) ? alpha*(-it->second) : 1.0;

			for (map::const_iterator it1 = A(k).begin(); it1 != A(k).end(); it1++)
			    B(i,it1->first) += v*it1->second;
		    }
		}
	    }
	}
	if (transform == TRANS_IU) {
	    /* A -> (I + alpha*U)*A */
	    alpha = 0.5;
	    for (uint i = 0; i < n; i++) {
		for (map::const_iterator it = A(i).begin(); it != A(i).end(); it++) {
		    uint k = it->first;
		    if (k >= i) {
			double v = (k != i) ? alpha*(-it->second) : 1.0;

			for (map::const_iterator it1 = A(k).begin(); it1 != A(k).end(); it1++)
			    B(i,it1->first) += v*it1->second;
		    }
		}
	    }
	}

	A_ = SkylineMatrix(B);
    }
}
