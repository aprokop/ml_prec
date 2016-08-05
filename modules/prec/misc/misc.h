#ifndef __MISC_H__
#define __MISC_H__

#include "include/exception.h"
#include "modules/matrix/matrix.h"
#include "modules/mesh/mesh.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

double cheb(double x, uint k);

void graph_planes(const std::string& filename, const SkylineMatrix& A, const uvector<uint>& map,
                  char plane, const MeshBase& mesh);

struct TailNode {
    uint index;
    /*
     * x[i] = a1*x[i+1] + F[i]
     * F[i] = a2*f[i] + a3*f[i-1]
     */
    double a1, a2, a3;
};
struct Tail : std::vector<TailNode> {
    char end_type;
    friend std::ostream& operator<<(std::ostream& os, const Tail& t) {
        os << t[0].index;
        for (uint i = 1; i < t.size(); i++)
            os << "-" << t[i].index;
        return os << "\n";
    }
};

enum LinkStatus {
    PRESENT,	/* Link is present */
    REMOVED,	/* Link is marked as removed */
    ABSENT	/* Link (i,j) was not in the original matrix */
};

class LinkTypeBase {
public:
    virtual ~LinkTypeBase() { }

    virtual void set_row(uint i) const = 0;

    /* Check link status */
    virtual LinkStatus stat(uint i, uint j) const = 0;
    virtual LinkStatus stat(uint j_) const = 0;
};

void construct_sparse_lu(const SkylineMatrix& A, const uvector<uint>& map, const uvector<uint>& rmap,
                         uint Md, uint M, const LinkTypeBase& ltype, double beta, const uvector<double>& aux,
                         SkylineMatrix& nA, SkylineMatrix& U, CSRMatrix& L);

#endif
