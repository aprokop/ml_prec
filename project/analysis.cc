#include <algorithm>
#include <numeric>
#include <fstream>

#include "main.h"
#include "include/logger.h"
#include "include/tools.h"

#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/solvers/solvers.h"

DEFINE_LOGGER("Analysis");

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

/*
 * Analyze the diagonal dominance in the matrix
 * For each row we calculate the ratio   (1-c/d)   , where c=sum(of all elements in row), and
 * d is the diagonal element.
 * Then, these ratios are used to construct the histogramm, using the precalculated baskets
 */
void anal_histogramm(const SkylineMatrix& A) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

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
    std::cout << "( ";
    for (uint i = 0; i+3 < n; i++)
	std::cout << bins[i] << ", ";
    std::cout << bins[n-3] << " )" << std::endl;
}

/*
 * Show number of dropped links for each value of q
 * Based on the current implementation of multi_split preconditioner (modules/multi_split)
 * For each value of q we perform procedure of the level 0 of the multi_split preconditioner
 * and calculate the number of removed links.
 * Additionally, we introduce parameter beta. It might be useful in the case when we have few
 * rows with very wead diagonal domination, and these rows completely distort the picture
 * The information is dumped in the file "qdropped.dat"
 */
void anal_qdropped(const SkylineMatrix& A) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    uvector<uint> sorted(1000);
    double c, d;
    double beta = 0.00;

    std::ofstream ofs("qdropped.dat");

    ofs << "0 0 " << std::endl;
    for (double q = 0.01; q < 0.99; q += 0.01) {
	std::cout << "q = " << q << std::endl;
	uint dropped = 0;

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

/*
 * Find value of q to remove a fixed number of links from each row
 * Based on the current implementation of multi_split preconditioner (modules/multi_split)
 * For each integer, 1...MAX_QS, we find the value of q, so that we are able to remove such
 * number of offdiagonal nonzeros per row.
 * Additionally, we introduce parameter beta. It might be useful in the case when we have few
 * rows with very wead diagonal domination, and these rows completely distort the picture
 */
void anal_q_rem_fixed_row(const SkylineMatrix& A) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    uint MAX_QS = 4; // try to remove up to MAX_QS number of links per row
    uvector<double> qs(MAX_QS, 0);

    std::ofstream ofs("iq.dat");

    uvector<uint> sorted(1000);
    double c, d;
    double beta = 0.1;
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

	double s = 0;
	uint k;
	for (k = 0; k < nrz && k < MAX_QS; k++) {
	    uint j_ = rstart+1 + sorted[k];

	    double aij = -a[j_];
	    s += aij;

	    double q = s/(c+s);

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
	std::cout << "To remove " << k+1 << " element(s) from each line we need q = " << qs[k] << std::endl;
}

/*
 * Analyze the magnitude of element difference inside rows
 * For each row we calculate the ratio  max_offdiagonal_element/min_diagonal_element
 * Then we take max and min for all rows
 */
void anal_offdiagonal_ratios(const SkylineMatrix& A) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

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
    }
    std::cout << "Ratio: [" << min << "," << max << "]\n";
}

/*
 * Mostly dumping procedure for raw values
 * For each row we calculate some quantity (for instance, 1-c/d, log10(c), c), and
 * then we dump it into "c_map.dat"
 */
void anal_unused1(const SkylineMatrix& A) {
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

    std::string filename("c_map.dat");
    std::ofstream os(filename.c_str(), std::ofstream::binary);
    os.write(reinterpret_cast<const char*>(&x[0]), N*sizeof(double));
}

/*
 * Creates and saves 1D z-matrices
 * This procedure assumes that the matrix A came from the standard 7-point
 * stencil discretization of the SPE test (it might be unsymmetric).
 * If we present matrix in the form A = D + A_xy + A_z, then the procedure
 * saves tridiagonal blocks of the matrix (D+A_z)
 * It is assumed, that the numeration is XYZ, i.e. index(i,j,k) = k*ny*nx + j*nx + k
 */
void anal_1D_Jacobi_array(const SkylineMatrix& A) {
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

		double s = A.row_sum(mesh_ind);

		v1 = v2 = 0.0;
		if (k != 0)	v1 = -a[ia[mesh_ind]+1];
		if (k != nz-1)	v2 = -a[ia[mesh_ind+1]-1];

		z_a[ind++] = s + (v1 + v2);
		if (k != 0)	z_a[ind++] = -v1;
		if (k != nz-1)	z_a[ind++] = -v2;
	    }
	    Az.set_size(nz,nz);

	    std::ostringstream oss;
	    oss << "z_matrices/matrix_i" << i << "_j" << j << ".mm";
	    dump(oss.str().c_str(), Az, MATRIX_MARKET);
	}
}

/*
 * Almost the same as anal_1D_Jacobi_array
 * The difference is that instead of the array we save the block-diagonal matrix,
 * with each block corresponding to a tridiagonal matrix
 * The resulting matrix is dumped into "az.mm"
 */
void anal_1D_Jacobi_global(const SkylineMatrix& A) {
    uint nx = 60, ny = 220, nz = 85;
    ASSERT(A.size() == 60*220*85, "Wrong matrix size");

    const uvector<uint> &ia = A.get_ia(), &ja = A.get_ja();
    const uvector<double> &a = A.get_a();

    SkylineMatrix Az;
    uvector<uint>& z_ia = Az.get_ia(), &z_ja = Az.get_ja();
    uvector<double>& z_a = Az.get_a();
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
	z_a[dind]  = A.row_sum(i);
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
}

/*
 * Calculate column dominance
 */
void anal_col_dominance(const SkylineMatrix& A) {
    uint N = A.size();
    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();
    const uvector<double>& a = A.get_a();

    const uint n = 8;
    double ticks[n] = { -1e6, -10, -1, 0, 1, 1.1, 10, 1e6};
    int bins[n] = {0, 0, 0, 0, 0, 0, 0, 0};

    uvector<double> x(N, 0.0);
    for (uint i = 0; i < N; i++)
	for (uint j = ia[i]; j < ia[i+1]; j++)
	    x[ja[j]] += a[j];
    for (uint i = 0; i < N; i++) {
	x[i] = 1 - x[i]/A(i,i);
	bins[std::upper_bound(ticks, ticks+n, x[i]) - (ticks+1)]++;
    }

    for (uint i = 0; i+1 < n; i++)
	std::cout << "[" << ticks[i] << "," << ticks[i+1] << ") : " << bins[i] << std::endl;
}

/*
 * Find the asymptotic convergence rate for nested iterations
 * We construct MultiSplitPrec with three levels, with different q for level 0 (q_0)
 * and level 1 (q_1). Level 2 is solved using direct method. We perform several
 * iterations with the preconditioner and study change in residual
 */
void anal_2level_convergence(const SkylineMatrix& A, const Config& cfg_) {
    Vector b(A.size()), x(A.size());

    const uint ntests = 3;
    const uint max_iter = 10;
    uvector<double> rates(ntests);

    uint   n = b.size();
    Vector r(n), z(n);
    double norm, init_norm;
    uint   niter;

    Config cfg = cfg_;
    cfg.sigmas.resize(2);
    cfg.max_levels = 3;
    cfg.unsym_matrix = true;
    cfg.prec = MULTI_SPLIT_PREC;
    cfg.niters = new_vector<uint>(1, 1);

    std::ofstream os("2level.dat");
    os << "# q = " << cfg.sigmas[0] << std::endl;

    srandom(time(NULL));

    for (double q1 = 0.0001; q1 < 0.99; q1 += 0.05) {
	cfg.sigmas[1] = q1;
	MultiSplitPrec B(A, cfg);
	LOG_DEBUG(B);

	for (uint i = 0; i < ntests; i++) {
	    LOG_INFO("q1 = " << q1 << ", TEST #" << i);

	    /* Generate x0 */
	    for (uint k = 0; k < x.size(); k++)
		x[k] = 20.*(random() - 0.5*RAND_MAX)/RAND_MAX + 100;

	    residual(A, b, x, r);

	    /* First iteration is atypical, ignore it */
	    B.solve(r, z);
	    daxpy(1., z, x);
	    residual(A, b, x, r);

	    norm = init_norm = calculate_norm(r, A, B, NORM_L2);

	    niter = 0;
	    while (niter < max_iter) {
		B.solve(r, z);
		daxpy(1., z, x);
		residual(A, b, x, r);

		norm = calculate_norm(r, A, B, NORM_L2);
		LOG_DEBUG("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

		niter++;
	    }

	    norm = calculate_norm(r, A, B, NORM_L2);

	    rates[i] = pow(norm/init_norm, 1./max_iter);
	}
	double rate = std::accumulate(rates.begin(), rates.end(), 0.0)/rates.size();

	std::cout << "q = " << cfg.sigmas[0] << ", q1 = " << q1 << ": rate = " << rate << std::endl;
	os << q1 << " " << rate << std::endl;
    }
    os.close();
}

void analyze(const SkylineMatrix& A, const Vector& b, const Config& cfg, AnalType analysis) {
    switch (analysis) {
	case ANAL_HISTOGRAMM	    : anal_histogramm(A); break;
	case ANAL_QDROPPED	    : anal_qdropped(A); break;
	case ANAL_Q_REM_FIXED_ROW   : anal_q_rem_fixed_row(A); break;
	case ANAL_OFFDIAGONAL_RATIOS: anal_offdiagonal_ratios(A); break;
	case ANAL_1D_JACOBI	    : anal_1D_Jacobi_global(A); break;
	case ANAL_COL_DOMINANCE	    : anal_col_dominance(A); break;
	case ANAL_2LEVEL_CONVERGENCE: anal_2level_convergence(A, cfg); break;
	case ANAL_NONE		    : LOG_WARN("Calling analysis with ANAL_NONE"); break;
    }
}

