#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "config/config.h"
#include "include/define.h"

#include <iostream>
#include <string>
#include <vector>

enum SolverType {
    PCG_SOLVER,
    CHEB_SOLVER,
    SIMPLE_SOLVER,
#ifdef HAVE_UMFPACK
    DIRECT_SOLVER
#endif
};

enum PrecType {
    UH_CHEB_PREC,	// Kuznetsov-Prokopenko
    AMG_PREC,		// AMG (Stüben)
    DIAG_PREC,		// Jacobi
    GS_PREC,		// Gauss-Seidel
#ifdef HAVE_UMFPACK
    BGS_PREC,		// block Gauss-Seidel
    RBGS_PREC,		// block Gauss-Seidel with a nonzero upper triangle ("relaxed" BGS)
#endif
    SYM_SPLIT_PREC,	// Kuznetsov-Prokopenko for symmetrized matrix
    MULTI_SPLIT_PREC	// Kuznetsov-Prokopenko (nested iterations for unsymmetric matrix)
};

enum AnalType {
    ANAL_NONE,
    ANAL_HISTOGRAMM,
    ANAL_QDROPPED,
    ANAL_Q_REM_FIXED_ROW,
    ANAL_OFFDIAGONAL_RATIOS,
    ANAL_1D_JACOBI
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

#endif // __CONFIG_H__
