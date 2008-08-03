#ifndef __MESH_H__
#define __MESH_H__

#include "modules/matrix/matrix.h"

struct Point {
    double x, y, z;
};

class Mesh {
private:
    static const uint	    nx = 60,	   ny = 220,	  nz = 85, N = nx*ny*nz;
    static const double	    hx = 20,	   hy = 10, 	  hz = 2;
    static const double size_x = 1200, size_y = 2200, size_z = 170;

    FVSparseMatrix A;
    std::vector<Point> nodes;

    uint index(uint i, uint j, uint k) const {
	return k*ny*nx + j*nx + i;
    }
    bool is_removed(uint i0, uint i1) const;

public:
    Mesh();

    void graph3D() const;
    void graph_xy_planes() const;
    void graph_z_lines() const;
};

#endif // __MESH_H__
