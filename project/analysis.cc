#include <algorithm>
#include <numeric>
#include <fstream>

#include "main.h"
#include "include/logger.h"
#include "include/tools.h"

DEFINE_LOGGER("Analysis");

void analysis(const SkylineMatrix& A) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    uvector<double> x(N);
    for (uint i = 0; i < N; i++) {
	double c = std::accumulate(a.begin() + ia[i], a.begin() + ia[i+1], 0.0);
	// x[i] = (log10(c) < -1) ? log10(c) : -1;
	x[i] = 1 - c/a[ia[i]];
    }
    // LOG_VAR(x);

    std::string filename("c_map.dat");
    std::ofstream os(filename.c_str(), std::ofstream::binary);
    os.write(reinterpret_cast<const char*>(&x[0]), N*sizeof(double));
}
