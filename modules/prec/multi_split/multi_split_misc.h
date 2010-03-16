#ifndef __MULTI_SPLIT_MISC_H__
#define __MULTI_SPLIT_MISC_H__

#include "include/define.h"
#include "include/exception.h"
#include "include/time.h"
#include "include/uvector.h"
#include "modules/matrix/matrix.h"

class LinkTypeMultiSplit {
private:
    uint n;
    const SkylineMatrix& A;
    uvector<char> a;

    /* Checks indices i and j (with ASSERT) and returns index in a */
    uint index(uint i, uint j) const THROW {
	ASSERT(i != j, "i == j == " << i);

	return A.index(i,j);
    }

public:
    LinkTypeMultiSplit(const SkylineMatrix& A_) : A(A_) {
	n = A_.size();
	a.resize(A.get_ja().size(), 1);
    }

    /* Mark directed link as removed from one end */
    bool mark(uint i, uint j) {
	return (--a[index(i,j)]) == 0;
    }
    /* Check whether the link was removed */
    bool is_removed(uint i, uint j) const {
	return a[index(i,j)] == 0;
    }
    /* Mark link as removed */
    void remove(uint i, uint j) {
	a[index(i,j)] = 0;
    }
};

#endif // __MULTI_SPLIT_MISC_H__
