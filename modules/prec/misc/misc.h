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

class LinkType {
private:
    std::vector<std::map<uint,char> > vec;

    void check(uint i, uint j) const {
	ASSERT(i != j, "i == j == " << i);
	ASSERT(i < vec.size() && j < vec.size(), "i = " << i << ", j = " << j << ", vec.size() == " << vec.size());
    }

public:
    LinkType(uint n) : vec(n) { }

    char operator()(uint i, uint j) const {
	check(i,j);

	uint i0, i1;
	if (i < j) { i0 = i; i1 = j; }
	else       { i0 = j; i1 = i; }

	std::map<uint,char>::const_iterator it = vec[i0].find(i1);
	if (it != vec[i0].end())
	    return it->second;
	return 0;
    }

    char& operator()(uint i, uint j) {
	check(i,j);
	if (i < j)
	    return vec[i][j];
	return vec[j][i];
    }
};

#endif
