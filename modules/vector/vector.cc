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
    std::ofstream os;
    if (type != BINARY) os.open(filename.c_str());
    else		os.open(filename.c_str(), std::ofstream::binary);

    const uint n = v.size();

    switch (type) {
	case ASCII : {
	    os << "# Vector size" << std::endl;
	    os << n << std::endl;
	    os << "# Vector values" << std::endl;
	    os << std::scientific << std::setprecision(15);
	    for (uint i = 0; i < n; i++)
		os << v[i] << std::endl;
	    break;
	}
	case HYPRE : {
	    os << "0 " << n-1 << std::endl;
	    os << std::scientific << std::setprecision(15);
	    for (uint i = 0; i < n; i++)
		os << i << " " << v[i] << std::endl;
	    break;
	}
        case BINARY: {
            os.write(reinterpret_cast<const char*>(&n), sizeof(uint));
            os.write(reinterpret_cast<const char*>(&v[0]), n*sizeof(uint));
            break;
        }
	default:
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
	default: THROW_EXCEPTION("Not implemented");
    }
    LOG_INFO("Loaded vector: size = " << n);
}

Vector vector_product(const Vector& v1, const Vector& v2) THROW {
    ASSERT_SIZE(v1.size(), 3);
    ASSERT_SIZE(v2.size(), 3);

    Vector v(3);
    v[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return v;
}

bool is_nan(const Vector& v) {
    uint n = v.size();
    for (uint i = 0; i < n; i++)
	if (::is_nan(v[i]))
	    return true;
    return false;
}
