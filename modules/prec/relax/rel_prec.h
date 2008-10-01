#ifndef __REL_PREC_H__
#define __REL_PREC_H__

#include "modules/prec/prec_base.h"
#include "modules/prec/misc/misc.h"
#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"
#include "modules/mesh/mesh.h"
#include <iostream>

class RelPrec : public PrecBase {
private:
    const Mesh& mesh;

    std::vector<double> sigmas;
    double gamma;
    uint niter;

    struct Level {
	uint N, nnz;

	mutable Vector x0, x1, f1;
	SkylineMatrix A;

	std::vector<uint> tr; 
	std::vector<uint> dtr;
	std::vector<Tail> tails;

	std::vector<double> aux;
    };
    uint nlevels;
    std::vector<Level> levels;

    void    construct_level(uint i, const SkylineMatrix& A);
    void    solve(Vector& f, Vector& x, uint level) THROW;

public:
    RelPrec(const SkylineMatrix& A, uint _niter, double _gamma, const std::vector<double>& _sigmas, const Mesh& _mesh);

    void graph_planes(const std::string& filename, uint level, char plane) const;
    void solve(Vector& f, Vector& x) THROW;
    friend std::ostream& operator<<(std::ostream& os, const RelPrec& p);
};

#endif
