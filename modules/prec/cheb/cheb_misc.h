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

enum LinkStatus {
    PRESENT,
    REMOVED,
    ABSENT
};

class LinkTypeCheb {
private:
    uint n;
    const SkylineMatrix& A;
    const uvector<uint>& ia;
    const uvector<uint>& ja;
    uvector<uint> a;

    uint row;

    /* Checks indices i and j (with ASSERT) and returns index in a */
    uint index(uint i, uint j) const THROW {
	ASSERT(i != j, "i == j == " << i);
	ASSERT(i < n && j < n, "i = " << i << ", j = " << j << ", n == " << n);

	if (i > j)
	    std::swap(i,j);

	return A.index(i,j);
    }

    uint real_j_(uint j_) const {
	ASSERT(ia[row] < j_ && j_ < ia[row+1], "Wrong j_: " << j_ << " (expected (" << ia[row] << ", " << ia[row+1] << ")");
	return (ja[j_] < row) ? a[j_] : j_;
    }


public:
    LinkTypeCheb(const SkylineMatrix& A_) : A(A_), ia(A_.get_ia()), ja(A_.get_ja()) {
	n = A.size();
	row = uint(-1);

	a.resize(ja.size());

	uint j_ = 0;
	for (uint i = 0; i < n; i++)
	    /* Set upper half to the link status and lower half to reference upper half */
	    for (j_ = ia[i]+1; j_ < ia[i+1]; j_++) {
		uint j = ja[j_];
		a[j_] = (j < i) ? index(j,i) : 2;
	    }
    }

    void set_row(uint i)    { row = i; }

    /* NOTE:
     * for mark(j_), stat(j_), remove(j_)
     * the general use case is the following:
     * if one wants to check/mark links in one row he should
     *	1. call set_row(row_index)
     *	2. access elements with the offset in A.ja array
     * Such procedure allows skip the search stage of index(i,j)
     */

    /* Mark link as removed from one end. Returns whether the link can be removed from both ends */
    LinkStatus mark(uint j_)		    {   return (--a[real_j_(j_)] == 0) ? REMOVED : PRESENT;   }
    LinkStatus mark(uint i, uint j) {
	uint ind = index(i,j);
	ASSERT(ind != uint(-1), "Cannot mark absent link: (" << i << "," << j << ")");
	return (--a[ind] == 0) ? REMOVED : PRESENT;
    }

    /* Check link status */
    LinkStatus stat(uint j_) const	    {   return a[real_j_(j_)] ? PRESENT : REMOVED;   }
    LinkStatus stat(uint i, uint j) const {
	uint ind = index(i,j);
	if (ind == uint(-1))
	    return ABSENT;
	return a[ind] ? PRESENT : REMOVED;
    }

    /* Mark link as removed */
    void remove(uint j_)		    {   a[real_j_(j_)] = 0;	  }
    void remove(uint i, uint j) {
	uint ind = index(i,j);
	ASSERT(ind != uint(-1), "Cannot removed absent link: (" << i << "," << j << ")");
	a[ind] = 0;
    }
};

#endif // __CHEB_MISC_H__
