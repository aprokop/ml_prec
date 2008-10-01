#ifndef __MESH_H__
#define __MESH_H__

#include "modules/matrix/matrix.h"

typedef unsigned int uint;

struct Point {
    double x, y, z;
};

#define ORDER_XYZ
// #define ORDER_ZXY
class Mesh {
public:
    const uint	 nx,	 ny,	 nz,	N;
    const double hx,	 hy, 	 hz;
    const double size_x, size_y, size_z;

    uint index(uint i, uint j, uint k) const {
	ASSERT(i < nx && j < ny && k < nz, "Wrong indices: (" << i << "," << j << "," << k << ")");

#if   defined ORDER_ZXY
	return j*nx*nz + i*nz + k; // z oriented ordering: yxz
#elif defined ORDER_XYZ
	return index_k(i, j, k);
#endif
    }


private:
    double c;

    std::vector<double> kx, ky, kz;
    std::vector<Point> nodes;

    uint index_k(uint i, uint j, uint k) const {
	return k*ny*nx + j*nx + i;
    }

public:
    Mesh(double _c);

    void construct_matrix(SkylineMatrix& A, uint nwells = 0) const;
    const std::vector<Point>& get_nodes() const {
	return nodes;
    }
};

#endif // __MESH_H__
