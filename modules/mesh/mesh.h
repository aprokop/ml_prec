#ifndef __MESH_H__
#define __MESH_H__

#include "include/define.h"
#include "modules/matrix/matrix.h"
#include <iostream>
#include <map>

struct Point {
    double x, y, z;
};

std::ostream& operator<<(std::ostream& os, const Point& p);

class MeshBase {
protected:
    uint   nx,	 ny,	 nz,	N;  /* Mesh dimensions */
    double hx,	 hy, 	 hz;	    /* Mesh steps */
    float  size_x, size_y, size_z;  /* Real mesh sizes */

    std::vector<Point> nodes;

public:
    virtual ~MeshBase() {}

    virtual uint index(uint i, uint j, uint k) const = 0;
    virtual void construct_matrix(SkylineMatrix& A, double c) const = 0;

    const std::vector<Point>& get_nodes() const {
        return nodes;
    }
    float get_size(char t) const {
        switch (t) {
            case 'x': return size_x;
            case 'y': return size_y;
            case 'z': return size_z;
            default: THROW_EXCEPTION("Unknown type '" << t << "'");
        }
    }
    uint get_n(char t) const {
        switch (t) {
            case 'x': return nx;
            case 'y': return ny;
            case 'z': return nz;
            default: THROW_EXCEPTION("Unknown type '" << t << "'");
        }
    }
};

#define XYZ
// #define ZXY
class SPEMesh : public MeshBase {
private:
    std::vector<double> kx, ky, kz;	    /* Diffusion coefficients */
    uint knx, kny, knz,	kN;

    uint index_k(uint i, uint j, uint k) const {
        return (k*kny*knx + j*knx + i) % kN;
    }

public:
    SPEMesh(uint _nx = 60, uint _ny = 220, uint _nz = 85);

    void construct_matrix(SkylineMatrix& A, double c) const;
    void construct_matrix_nonsym(SkylineMatrix& A, double c, double shift) const;

    uint index(uint i, uint j, uint k) const{
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
