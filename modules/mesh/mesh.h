#ifndef __MESH_H__
#define __MESH_H__

#include "modules/matrix/matrix.h"

typedef unsigned int uint;

struct Point {
    double x, y, z;
};

class Mesh {
private:
    const uint	 nx,	 ny,	 nz,	N;
    const double hx,	 hy, 	 hz;
    const double size_x, size_y, size_z;

    double c;
    std::vector<double> kx, ky, kz;
    std::vector<Point> nodes;

    uint index_k(uint i, uint j, uint k) const {
	return k*ny*nx + j*nx + i;
    }

    // x oriented ordering: zyx
    uint index(uint i, uint j, uint k) const {
	ASSERT(i < nx && j < ny && k < nz, "Wrong indices: (" << i << "," << j << "," << k << ")");
#if 1
	return j*nx*nz + i*nz + k; // z oriented ordering: yxz
#else
	return index_k(i, j, k);
#endif
    }

public:
    Mesh(double _c);

    void construct_matrix(SkylineMatrix& A, uint nwells = 0) const;

    friend class Prec;
};

#endif // __MESH_H__
