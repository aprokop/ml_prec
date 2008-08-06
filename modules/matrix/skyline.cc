#include "matrix.h"
#include "include/logger.h"
#include "include/exception.h"
#include "include/tools.h"

DEFINE_LOGGER("SkylineMatrix");

SkylineMatrix::SkylineMatrix() {
    nrow = ncol = 0;
}

double SkylineMatrix::get(uint i, uint j) const THROW {
    check_indices(i, j);

    // in case matrix is in skyline mode
    if (i == j)
	return a[ia[i]];
    std::vector<uint>::const_iterator start = ja.begin() + ia[i] + 1; 
    std::vector<uint>::const_iterator   end = ja.begin() + ia[i+1];
    std::vector<uint>::const_iterator it = std::lower_bound(start, end, j);
    if (it != end && !(*it < j)) 
	return a[it-ja.begin()];

    LOG_WARN("Returning zero element for i = " << i << ", j = " << j);
    return 0;
}

void SkylineMatrix::add(uint i, uint j, double v) THROW {
    check_indices(i, j);

    // in case matrix is in skyline form
    if (i == j) {
	a[ia[i]] += v;
	return;
    }

    std::vector<uint>::iterator start = ja.begin() + ia[i] + 1; 
    std::vector<uint>::iterator   end = ja.begin() + ia[i+1];
    std::vector<uint>::iterator it = std::lower_bound(start, end, j);
    if (it == end || j < *it) {
	// creating new element
	uint pos = it - ja.begin();
	ja.insert(it, j);
	a.insert (a.begin() + pos, v);
	for (uint k = i+1; k <= nrow; k++)
	    ia[k]++;
    } else {
	// adding to existing element
	a[it - ja.begin()] += v;
    }
}
