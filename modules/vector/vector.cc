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

void dump(const Vector& v, const std::string& filename, bool ascii) {
    if (ascii == false) {
	THROW_EXCEPTION("Not implemented");
    } else {
	std::ofstream os(filename.c_str());
	const uint n = v.size();

	os << "# Vector size" << std::endl;
	os << n << std::endl;
	os << "# Vector values" << std::endl;
	os << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::scientific;
	for (uint i = 0; i < n; i++)
	    os << v[i] << std::endl;
    }
}

void load(Vector& v, const std::string& filename, bool ascii) {
    uint n;
    if (ascii == false) {
	THROW_EXCEPTION("Not implemented");
    } else {
	std::ifstream is(filename.c_str());
	ASSERT(is.good(), "Problem reading file \"" << filename  << "\"");

	const int SN = 2009;
	char str[SN];

	is.getline(str, SN); /* "# Vector size" */
	is >> n;
	v.resize(n);

	is.getline(str, SN);
	is.getline(str, SN); /* "# Vector values" */
	for (uint i = 0; i < n; i++)
	    is >> v[i];
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
