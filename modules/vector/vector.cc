#include "vector.h"
#include "include/logger.h"
#include "include/tools.h"

#include "config/config.h"
#ifdef HAVE_BLAS
#include "include/blas.h"
#endif

#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

DEFINE_LOGGER("Vector");

void dump(const std::string& filename, const Vector& v, DumpType type) THROW {
    switch (type) {
	case ASCII : {
	    std::ofstream os(filename.c_str());
	    const uint n = v.size();

	    os << "# Vector size" << std::endl;
	    os << n << std::endl;
	    os << "# Vector values" << std::endl;
	    os << std::scientific << std::setprecision(15);
	    for (uint i = 0; i < n; i++)
		os << v[i] << std::endl;
	    break;
	}
	case HYPRE : {
	    std::ofstream os(filename.c_str());
	    const uint n = v.size();

	    os << "0 " << n-1 << std::endl;
	    os << std::scientific << std::setprecision(15);
	    for (uint i = 0; i < n; i++)
		os << i << " " << v[i] << std::endl;
	    break;
	}
	case BINARY:
	    THROW_EXCEPTION("Not implemented");
    }
}

void load(Vector& v, const std::string& filename, DumpType type) THROW {
    uint n;
    switch (type) {
	case ASCII : {
	    std::ifstream is(filename.c_str());
	    if (!is.good())
		THROW_EXCEPTION("Problem reading file \"" << filename  << "\"");

	    const int SN = 2009;
	    char str[SN];

	    is.getline(str, SN); /* "# Vector size" */
	    is >> n;
	    v.resize(n);

	    is.getline(str, SN);
	    is.getline(str, SN); /* "# Vector values" */
	    for (uint i = 0; i < n; i++)
		is >> v[i];
	    break;
	}
	case HYPRE :
	case BINARY:
		    THROW_EXCEPTION("Not implemented");
    }
    LOG_INFO("Loaded vector: size = " << n);
}

bool is_nan(const Vector& v) {
    uint n = v.size();
    for (uint i = 0; i < n; i++)
	if (::is_nan(v[i]))
	    return true;
    return false;
}
