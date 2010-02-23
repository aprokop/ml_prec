#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "include/define.h"

#include <iostream>
#include <string>
#include <vector>

enum SolverType {
    PCG_SOLVER,
    CHEB_SOLVER
};

enum PrecType {
    UH_CHEB_PREC,
    AMG_PREC,
    DIAG_PREC
};

struct Config {
    /* Chebyshev iterations parameters */
    std::vector<uint>   niters;
    std::vector<double> sigmas;

    /* Mesh parameters */
    double c;
    uint   nx, ny, nz;

    /* Run parameters */
    uint   ntests;
    bool   use_tails;
    bool   optimize_storage;
    bool   unsym_matrix;
    SolverType solver;
    PrecType prec;

    std::string matrix;

    friend std::ostream& operator<<(std::ostream& os, const Config& config);
};

double cheb(double x, uint k);

#endif // __CONFIG_H__
