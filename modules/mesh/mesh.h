#ifndef __MESH_H__
#define __MESH_H__

#include "include/define.h"
#include "modules/matrix/matrix.h"
#include <map>

struct Point {
    double x, y, z;
};

class MeshBase {
protected:
    uint   nx,	 ny,	 nz,	N;
    double hx,	 hy, 	 hz;
    float  size_x, size_y, size_z; /* Needed for graph_planes */

    std::vector<Point> nodes;

public:
    virtual ~MeshBase() {}

    virtual uint index(uint i, uint j, uint k) const = 0;
    virtual void construct_matrix(SkylineMatrix& A, double c) const = 0;

    friend void graph_planes(const std::string& filename, const SkylineMatrix& A,
			     const std::map<uint,uint>& rev_map,
			     char plane, bool map_identity, const MeshBase& mesh);
};

#define XYZ
// #define ZXY
class SPEMesh : public MeshBase {
private:
    std::vector<double> kx, ky, kz;
    uint knx, kny, knz,	kN;

    uint index_k(uint i, uint j, uint k) const {
	return (k*kny*knx + j*knx + i) % kN;
    }

public:
    SPEMesh(uint _nx = 60, uint _ny = 220, uint _nz = 85);

    void construct_matrix(SkylineMatrix& A, double c) const;
    void construct_matrix_unsym(SkylineMatrix& A, double c, double shift) const;

    uint index(uint i, uint j, uint k) const {
	ASSERT(i < nx && j < ny && k < nz, "Wrong indices: (" << i << "," << j << "," << k << ")");

#ifdef XYZ
	/* XYZ */
	return k*ny*nx + j*nx + i;
#elif defined ZXY
	/* ZXY */
	return j*nx*nz + i*nz + k;
#endif
    }
};

#endif // __MESH_H__
