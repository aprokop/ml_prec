#include "matrix.h"
#include "include/logger.h"

#include <iostream>

DEFINE_LOGGER("FVSparseMatrix");

void FVSparseMatrix::new_link(uint i0, uint i1, double x) THROW {
    // LOG_DEBUG("Adding new link: (" << i0 << "," << i1 << ") : " << x);

    ASSERT(i0 < nrow && i0 < ncol && i1 < nrow && i1 < ncol && i0 != i1, 
	   "Wrong indices: i0 = " << i0 << ", i1 = " << i1 << ", nrow = " << nrow << ", ncol = " << ncol);
    ASSERT(x > 0, "Entry must be positive");
    ASSERT(vrows[i0].find(i1) == vrows[i0].end(), "Link is already present.");

    vrows[i0][i0] +=  x;
    vrows[i0][i1]  = -x;
    vrows[i1][i0]  = -x;
    vrows[i1][i1] +=  x;
}

double FVSparseMatrix::remove_link(uint i0, uint i1) THROW {
    // LOG_DEBUG("Removing link: (" << i0 << "," << i1 << ")");

    ASSERT(i0 < nrow && i0 < ncol && i1 < nrow && i1 < ncol && i0 != i1, 
	   "Wrong indices: i0 = " << i0 << ", i1 = " << i1 << ", nrow = " << nrow << ", ncol = " << ncol);

    Row::iterator it = vrows[i0].find(i1);
    if (it == vrows[i0].end())
	THROW_EXCEPTION("There is no link [" << i0 << ", " << i1 << "] in matrix");

    double d = it->second; // d is negative
    vrows[i0].erase(i1);
    vrows[i1].erase(i0);
    vrows[i0][i0] += d;
    vrows[i1][i1] += d;

    return -d;
}

