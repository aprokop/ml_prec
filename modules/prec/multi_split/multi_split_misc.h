#ifndef __MULTI_SPLIT_MISC_H__
#define __MULTI_SPLIT_MISC_H__

#include "include/define.h"
#include "include/exception.h"
#include "include/time.h"
#include "include/uvector.h"
#include "modules/matrix/matrix.h"

class LinkTypeMultiSplit : public LinkTypeBase {
private:
    uint n;
    const SkylineMatrix& A;
    uvector<char> a;

    // checks indices i and j (with ASSERT) and returns index in a
    uint index(uint i, uint j) const{
        ASSERT(i != j, "i == j == " << i);
        return A.index(i,j);
    }

public:
    LinkTypeMultiSplit(const SkylineMatrix& A_) : A(A_) {
        n = A.size();
        a.resize(A.get_ja().size(), 1);
    }

    // check link status
    LinkStatus stat(uint _j) const {   return a[_j] ? PRESENT : REMOVED;   }
    LinkStatus stat(uint i, uint j) const {
        uint ind = index(i,j);
        if (ind == uint(-1))
            return ABSENT;
        return stat(ind);
    }

    // extra function not used here
    void set_row(uint i) const { }

    // mark link as removed
    void remove(uint _j) {   a[_j] = 0;   }
    void remove(uint i, uint j) {
        uint ind = index(i,j);
        ASSERT(ind != uint(-1), "Cannot removed absent link: (" << i << "," << j << ")");
        a[ind] = 0;
    }

    // Find the remaining link from a point (actually it finds the first one).
    // Returns index in ja
    // The caller must be sure that only one link remains
    bool find_remaining_link(uint i, uint& _j) {
        ASSERT(i < n, "i = " << i << ", n = " << n);

        const uvector<uint>& ia = A.get_ia();
        const uvector<uint>& ja = A.get_ja();

        for (_j = ia[i]+1; _j < ia[i+1]; _j++)
            if (a[_j])
                return true;

        return false;
    }

    // Remove remaining link from a point
    // As in find_remaining_link, it actually removes the first one
    // The caller must be sure that only one link remains
    bool remove_remaining_link(uint i, uint& _j) {
        bool ret = find_remaining_link(i, _j);
        if (ret)
            a[_j] = 0;
        return ret;
    }
};

#endif // __MULTI_SPLIT_MISC_H__
