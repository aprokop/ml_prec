#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "include/define.h"

#include <iostream>
#include <string>
#include <vector>

enum SolverType {
    PCG_SOLVER,
    CHEB_SOLVER,
    SIMPLE_SOLVER,
    DIRECT_SOLVER
};

enum PrecType {
    UH_CHEB_PREC,	// Kuznetsov-Prokopenko
    AMG_PREC,		// AMG (Stüben)
    DIAG_PREC,		// Jacobi
    GS_PREC,		// Gauss-Seidel
    SYM_SPLIT_PREC,	// Kuznetsov-Prokopenko for symmetrized matrix
    MULTI_SPLIT_PREC	// Kuznetsov-Prokopenko (nested iterations for unsymmetric matrix)
};

enum AnalType {
    ANAL_NONE,
    ANAL_HISTOGRAMM,
    ANAL_QDROPPED,
    ANAL_Q_REM_FIXED_ROW
};

struct Config {
    /* Chebyshev iterations parameters */
    std::vector<uint>   niters;		/* Number of level iterations */
    std::vector<double> sigmas;		/* Level sigmas (or q for unsymmetric) */

    /* Mesh parameters */
    double c;				/* Reaction coefficient */
    uint   nx, ny, nz;			/* Mesh dimensions */

    /* Run parameters */
    uint   ntests;			/* Number of tests for averaging */
    bool   use_tails;			/* Optional use of tail elimination */
    bool   optimize_storage;		/* Storage optimization for symmetric matrices */
    bool   unsym_matrix;		/* Symmetricity of the matrix */
    double unsym_shift;			/* Degree of the generated unsymmetry */
    SolverType solver;			/* External solver */
    PrecType prec;			/* Preconditioner */

    /* Other */
    bool   dump_data;			/* Dumping matrix and rhs */
    std::string dir;			/* Directory to dump the results */
    AnalType analysis;			/* Type of matrix analysis */

    std::string matrix;			/* Mesh file */
    std::string vector;			/* Vector file */

    friend std::ostream& operator<<(std::ostream& os, const Config& config);
};

double cheb(double x, uint k);

#endif // __CONFIG_H__
