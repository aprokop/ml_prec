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
    const MeshBase& mesh;

    std::vector<double> sigmas;
    double gamma;
    uint niter;

    struct Level {
        uint N, nnz;

        mutable Vector x1, f1;
        Vector tmp0, tmp1;
        SkylineMatrix A;

        std::vector<uint> tr;
        std::vector<uint> dtr;
        std::vector<Tail> tails;

        Vector aux;
    };
    uint nlevels;
    std::vector<Level> levels;

    void    construct_level(uint i, const SkylineMatrix& A);
    void    solve(Vector& f, Vector& x, uint level);

public:
    RelPrec(const SkylineMatrix& A, uint _niter, double _gamma, const std::vector<double>& _sigmas, const MeshBase& _mesh);

    void solve(Vector& f, Vector& x);
    void graph_planes(const std::string& filename, uint level, char plane) const;

    friend std::ostream& operator<<(std::ostream& os, const RelPrec& p);
};

#endif
