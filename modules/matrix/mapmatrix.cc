#include "matrix.h"
#include "include/logger.h"
#include "include/tools.h"

DEFINE_LOGGER("MapMatrix");

std::ostream& operator<<(std::ostream& os, const MapMatrix& sm) {
    os << "Size: " << sm.nrow << "x" << sm.ncol << std::endl;
    for (uint i = 0; i < sm.nrow; i++) {
	os << "  Row: " << i << std::endl;
	const MapMatrix::Row& row = sm.data[i];
	for (MapMatrix::Row::const_iterator it = row.begin(); it != row.end(); it++)
	    os << "    " << it->first << ": " << it->second << std::endl;
	os << std::endl;
    }

    return os;
}
