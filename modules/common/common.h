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
    std::string matrix;
    SolverType solver;
    
    friend std::ostream& operator<<(std::ostream& os, const Config& config);
};

double cheb(double x, uint k);

#endif // __CONFIG_H__
