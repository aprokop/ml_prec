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

void analysis_1D_Jacobi(const SkylineMatrix& A) {
    uint nx = 60, ny = 220, nz = 85;
    ASSERT(A.size() == 60*220*85, "Wrong matrix size");

    const uvector<uint> &ia = A.get_ia(), &ja = A.get_ja();
    const uvector<double> &a = A.get_a();

    SkylineMatrix Az;
    uvector<uint>& z_ia = Az.get_ia(), &z_ja = Az.get_ja();
    uvector<double>& z_a = Az.get_a();

    /* Construct tridiagonal stencil */
    z_ia.resize(nz+1);
    z_ja.resize(3*nz);
    z_ia[0] = 0;
    uint ind = 0;
    for (uint k = 0; k < nz; k++) {
	z_ja[ind++] = k;
	if (k != 0)	z_ja[ind++] = k-1;
	if (k != nz-1)	z_ja[ind++] = k+1;

	z_ia[k+1] = ind;
    }
    z_ja.resize(ind);

    /* Fill out values */
    z_a.resize(ind);
    for (uint j = 0; j < ny; j++)
	for (uint i = 0; i < nx; i++) {
	    double v1, v2;
	    ind = 0;
	    for (uint k = 0; k < nz; k++) {
		uint mesh_ind = k*ny*nx + j*nx + i;

		double s = 0.0;
		for (uint l = ia[mesh_ind]; l < ia[mesh_ind+1]; l++)
		    s += a[l];

		v1 = v2 = 0.0;
		if (k != 0)	v1 = -a[ia[mesh_ind]+1];
		if (k != nz-1)	v2 = -a[ia[mesh_ind+1]-1];

		z_a[ind++] = s + (v1 + v2);
		if (k != 0)	z_a[ind++] = -v1;
		if (k != nz-1)	z_a[ind++] = -v2;
	    }
	    Az.set_size(nz,nz);

	    {
		/* Dump matrix in MatrixMarket ASCII format */
		std::ostringstream oss;
		oss << "z_matrices/matrix_i" << i << "_j" << j << ".mm";
		std::ofstream os(oss.str().c_str());
		os << std::fixed << std::setprecision(15);

		os << "%%MatrixMarket matrix coordinate real general" << std::endl;
		os << nz << " " << nz << " " << ind << std::endl;
		for (uint p = 0; p < nz; p++)
		    for (uint q = z_ia[p]; q < z_ia[p+1]; q++)
			os << p+1 << " " << z_ja[q]+1 << " " << z_a[q] << std::endl;
	    }

#if 0
	    /* Check matrix */
	    for (uint k = 1; k < nz-1; k++) {
		uint i0 = k*ny*nx + j*nx + i;
		uint i1 = (k-1)*ny*nx + j*nx + i;
		uint i2 = (k+1)*ny*nx + j*nx + i;
		if (!is_equal(A(i0,i1), Az(k,k-1)) ||
		    !is_equal(A(i0,i2), Az(k,k+1)))
		    THROW_EXCEPTION("Smth is wrong at k = " << k << ", (" << i << "," << j << ")");
	    }
#endif
	}
}
