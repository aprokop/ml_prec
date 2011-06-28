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

#if 0
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

	}
#else
    /* Construct Az matrix */
    uint n = A.size();
    z_ia.resize(n+1);
    z_ja.resize(3*n);
    z_a.resize(3*n);

    double ind = 0, dind = 0;
    z_ia[0] = 0;
    for (uint i = 0; i < n; i++) {
	dind = ind;
	z_ja[dind] = i;
	z_a[dind]  = std::accumulate(&a[ia[i]], &a[ia[i+1]], 0.0);
	ind++;
	for (uint j = ia[i]+1; j < ia[i+1]; j++)
	    if (((abs(ja[j] - i)) % (nx*ny)) == 0) {
		z_ja[ind]  = ja[j];
		z_a[ind]   = a[j];
		z_a[dind] -= a[j];
		ind++;
	    }

	z_ia[i+1] = ind;
    }
    z_ja.resize(ind);
    z_a.resize(ind);
    Az.set_size(n, n);

    dump("az.mm", Az, MATRIX_MARKET);
#endif
}

/*
 * Insertion sort of vector sorted with respect to absolute values in vector a
 * NOTE: later if for some matrices we would get many elements in a row we could
 * replace this sort with a faster one (think heapsort)
 * T must support fabs and < operators
 */
static void psort(const double *a, uint n, uvector<uint>& sorted) {
    sorted[0] = 0;
    double v;
    int j;
    for (uint i = 1; i < n; i++) {
	v = fabs(a[i]);
	for (j = i-1; j >= 0 && fabs(a[sorted[j]]) > v; j--)
	    sorted[j+1] = sorted[j];
	sorted[j+1] = i;
    }
}

void examine(const SkylineMatrix& A, AnalType analysis) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    if (analysis == ANAL_Q_REM_FIXED_ROW) {
	/* Find qs to remove a fixed number of links from each row */
	uint MAX_QS = 4;
	uvector<double> qs(MAX_QS, 0);

	uint IQ = 1;
	std::ofstream ofs("iq.dat");

	uvector<uint> sorted(1000);
	double c, d;
	double beta = 0.1;
	for (uint i = 0; i < N; i++) {
	    uint rstart = ia[i];		/* Row start */
	    uint rend   = ia[i+1];		/* Row end */
	    uint nrz    = rend - rstart - 1;    /* Number of outgoing links */

	    d = a[rstart]; // diagonal element
	    c = 0;
	    for (uint j_ = rstart; j_ < rend; j_++)
		c += a[j_];

	    if (1-c/d > 1-beta)
		continue;

	    /* Sort off-diagonal elements wrt their abs values */
	    psort(&a[0] + rstart+1, nrz, sorted);

	    double s = 0;
	    uint k;
	    for (k = 0; k < nrz && k < MAX_QS; k++) {
		uint j_ = rstart+1 + sorted[k];

		double aij = -a[j_];
		s += aij;

		double q = s/(c+s);

		// if (k == IQ-1 && q > 0.0001)
		    // ofs << i << " " << q << std::endl;

		if (q > qs[k])
		    qs[k] = q;
	    }
	    if (k == nrz)
		for (; k < MAX_QS; k++) {
		    double q = s/(c+s);
		    if (q > qs[k])
			qs[k] = q;
		}
	}
	std::cout << std::fixed << std::setprecision(3);
	for (uint k = 0; k < MAX_QS; k++)
	    std::cout << "To remove " << k+1 << " element from each line we need q = " << qs[k] << std::endl;
    }
    if (analysis == ANAL_QDROPPED) {
	/* Show number of dropped links for each value of q */
	uvector<uint> sorted(1000);
	double c, d;
	double beta = 0.00;

	std::ofstream ofs("qdropped.dat");

	ofs << "0 0 " << std::endl;
	for (double q = 0.01; q < 0.99; q += 0.01) {
	    std::cout << "q = " << q << std::endl;
	    uint dropped = 0;

	    for (uint i = 0; i < N; i++) {
		uint rstart = ia[i];		    /* Row start */
		uint rend   = ia[i+1];		    /* Row end */
		uint nrz    = rend - rstart - 1;    /* Number of outgoing links */

		d = a[rstart]; // diagonal element
		c = 0;
		for (uint j_ = rstart; j_ < rend; j_++)
		    c += a[j_];

		if (1-c/d > 1-beta)
		    continue;

		/* Sort off-diagonal elements wrt their abs values */
		psort(&a[0] + rstart+1, nrz, sorted);

		double s = q/(1-q)*c;
		for (uint k = 0; k < nrz; k++) {
		    uint j_ = rstart+1 + sorted[k];

		    double aij = -a[j_];
		    if (aij <= s) {
			s -= aij;
			dropped++;
		    }
		}
	    }
	    ofs << q << " " << dropped << std::endl;
	}
	ofs << "1.0 " << ia[N]-N << std::endl;
    }
    if (0) {
	uvector<uint> sorted(1000);
	double min = 1e50, max = 0;
	for (uint i = 0; i < N; i++) {
	    uint rstart = ia[i];		    /* Row start */
	    uint rend   = ia[i+1];		    /* Row end */
	    uint nrz    = rend - rstart - 1;    /* Number of outgoing links */

	    /* Sort off-diagonal elements wrt their abs values */
	    psort(&a[0] + rstart+1, nrz, sorted);

	    double ratio = a[rstart+1+sorted[nrz-1]]/a[rstart+1+sorted[0]];
	    if (ratio < min) min = ratio;
	    if (ratio > max) max = ratio;
	    std::cout << i << ": " << ratio<< std::endl;
	}
	std::cout << "Ratio: [" << min << "," << max << "]\n";
    }
    if (analysis == ANAL_HISTOGRAMM) {
	const uint n = 14;
	double ticks[n] = { 0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.9, 1, 1.0000001 };
	int bins[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	double sum, d;
	for (uint i = 0; i < N; i++) {
	    sum = 0;
	    d = a[ia[i]];
	    for (uint j = ia[i]; j < ia[i+1]; j++)
		sum += a[j];
	    bins[std::upper_bound(ticks, ticks+n, 1-sum/d) - (ticks+1)]++;
	}
	for (uint i = 0; i+1 < n; i++)
	    std::cout << "[" << ticks[i] << "," << ticks[i+1] << ") : " << bins[i] << std::endl;
    }
}
