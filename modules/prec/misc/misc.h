#ifndef __MISC_H__
#define __MISC_H__

#include "include/exception.h"
#include "include/tools.h"
#include "modules/matrix/matrix.h"
#include "modules/mesh/mesh.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

void graph_planes(const std::string& filename, const SkylineMatrix& A, const std::map<uint,uint>& rev_map,
		  char plane, bool map_identity, const MeshBase& mesh);

struct TailNode {
    uint index;
    // x_i = a1*x_{i+1} + a2*f_i + a3*f_{i-1}
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

#endif
