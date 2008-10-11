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
    double size_x, size_y, size_z;

    std::vector<Point> nodes;

public:
    virtual ~MeshBase() {}

    virtual uint index(uint i, uint j, uint k) const = 0; 
    virtual void construct_matrix(SkylineMatrix& A, double c) const = 0;

    friend void graph_planes(const std::string& filename, const SkylineMatrix& A, 
			     const std::map<uint,uint>& rev_map,
			     char plane, bool map_identity, const MeshBase& mesh);
};

#define ORDER_XYZ
// #define ORDER_ZXY
class SPEMesh : public MeshBase {
private:
    std::vector<double> kx, ky, kz;

    uint index_k(uint i, uint j, uint k) const {
	return k*ny*nx + j*nx + i;
    }

public:
    SPEMesh();

    void construct_matrix(SkylineMatrix& A, double c) const;
    
    uint index(uint i, uint j, uint k) const {
	ASSERT(i < nx && j < ny && k < nz, "Wrong indices: (" << i << "," << j << "," << k << ")");

#if   defined ORDER_ZXY
	return j*nx*nz + i*nz + k; // z oriented ordering: yxz
#elif defined ORDER_XYZ
	return index_k(i, j, k);
#endif
    }
};

#endif // __MESH_H__
