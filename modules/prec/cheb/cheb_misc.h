#ifndef __CHEB_MISC_H__
#define __CHEB_MISC_H__

#include "include/define.h"
#include "include/exception.h"
#include "include/time.h"
#include "include/uvector.h"
#include "modules/matrix/matrix.h"

#include <algorithm>
#include <map>
#include <vector>

class LinkType {
private:
    uint n;
    const SkylineMatrix& A;
    uvector<char> a;

    /* Checks indices i and j (with ASSERT) and returns index in a */
    uint index(uint i, uint j) const THROW {
	ASSERT(i != j, "i == j == " << i);
	ASSERT(i < n && j < n, "i = " << i << ", j = " << j << ", n == " << n);

	if (i > j)
	    std::swap(i,j);

	return A.index(i,j);
    }

public:
    LinkType(const SkylineMatrix& _A) : A(_A) {
	n = _A.size();
	a.resize(A.get_ja().size(), 2);
    }

    /* Mark link as removed from one end. Returns whether the link can be removed from both ends */
    bool mark(uint i, uint j) {
	return --a[index(i,j)] == 0;
    }
    /* Check whether the link was removed */
    bool is_removed(uint i, uint j) const {
	return a[index(i,j)] == 0;
    }
    /* Mark link as removed */
    void remove(uint i, uint j) {
	a[index(i,j)] = 0;
    }
    /* Find the remaining link (actually it finds the first one). Returns index in ja */
    bool remove_remaining_link(uint i, uint& _j) {
	ASSERT(i < n, "i = " << i << ", n = " << n);

	const uvector<uint>& ia = A.get_ia();
	const uvector<uint>& ja = A.get_ja();

	for (_j = ia[i+1] - 1; _j >= ia[i]+1; _j--) {
	    uint j = ja[_j];
	    if (j > i && a[_j]) {
		a[_j] = 0;
		return true;
	    }
	    if (j < i) {
		uint ind = index(i, j);
		if (a[ind]) {
		    a[ind] = 0;
		    return true;
		}
	    }
	}
	return false;
    }
};

#endif // __CHEB_MISC_H__
