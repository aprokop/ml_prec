#ifndef __PREC_H__
#define __PREC_H__

#include "modules/prec/misc/misc.h"
#include "modules/prec/prec_base.h"
#include "modules/mesh/mesh.h"
#include <iostream>

class Prec : public PrecBase {
private:
    const MeshBase& mesh;

    double galpha, gbeta;

    struct Level {
	uint N, nnz;

	mutable Vector x1, f1, u0, u1;
	SkylineMatrix A; // is not set for level 0

	std::vector<uint> tr; 
	std::vector<uint> dtr;
	std::vector<Tail> tails;

	std::vector<double> aux;

	double lmin, lmax;
    };
    uint nlevels;
    std::vector<Level> levels;
    uint ncheb;

    void    construct_level(uint i, const SkylineMatrix& A);
    void    solve(Vector& f, Vector& x, uint level) THROW;

public:
    Prec(double eps, uint ncheb, const SkylineMatrix& A, const MeshBase& _mesh);

    void solve(Vector& f, Vector& x) THROW;
    void graph_planes(const std::string& filename, uint level, char plane) const;

    double lmin() const { return levels[0].lmin; }
    double lmax() const { return levels[0].lmax; }

    friend std::ostream& operator<<(std::ostream& os, const Prec& p);
};

#endif
