#include "include/logger.h"
#include "include/time.h"
#include "include/tools.h"
#include "mesh.h"

#include <fstream>
#include <algorithm>

DEFINE_LOGGER("SPEMesh");

SPEMesh::SPEMesh(uint _nx, uint _ny, uint _nz) {
    ASSERT(_nx && _ny && _nz, "Wrong mesh dimensions: " << _nx << " x " << _ny << " x " << _nz);
    knx = 60;	    kny = 220;	    knz = 85;	    kN = knx*kny*knz;
    hx  = 20;	    hy  = 10;	    hz  = 2;
    nx  = _nx;      ny  = _ny;      nz  = _nz;

    N = nx*ny*nz;

    size_x = nx*hx;
    size_y = ny*hy;
    size_z = nz*hz;

    kx.resize(kN);
    ky.resize(kN);
    kz.resize(kN);


#if 1
    std::vector<double> kxyz;
    kxyz.reserve(3*kN);

    std::ifstream is("spe_perm.dat");

    if (!is.good())
        THROW_EXCEPTION("Problem reading file \"spe_perm.dat\"");

    std::copy(std::istream_iterator<double>(is), std::istream_iterator<double>(),
              std::back_inserter(kxyz));
    ASSERT(kxyz.size() == 3*kN, "Wrong spe_perm.dat size");

    memcpy(&kx[0], &kxyz[   0], kN*sizeof(double));
    memcpy(&ky[0], &kxyz[  kN], kN*sizeof(double));
    memcpy(&kz[0], &kxyz[2*kN], kN*sizeof(double));
#else
    std::ifstream is("spe_perm.dat", std::ifstream::binary);

    if (!is.good())
        THROW_EXCEPTION("Problem reading file \"spe_perm.dat\"");

    // Binary read (the file must have been converted first)
    is.read(reinterpret_cast<char*>(&kx[0]),  kN*sizeof(double));
    is.read(reinterpret_cast<char*>(&ky[0]),  kN*sizeof(double));
    is.read(reinterpret_cast<char*>(&kz[0]),  kN*sizeof(double));
#endif

    /* Construct nodes */
    nodes.resize(N);
    for (uint k = 0; k < nz; k++)
        for (uint j = 0; j < ny; j++)
            for (uint i = 0; i < nx; i++) {
                uint ind = index(i,j,k);
                nodes[ind].x = i*hx;
                nodes[ind].y = j*hy;
                nodes[ind].z = k*hz;
            }
    LEAVE_MESSAGE("Nodes constructed");
}

void SPEMesh::construct_matrix(SkylineMatrix& A, double c) const {
    LOG_INFO("Constructing symmetric part");

    uint i0, i1;
    double v;

    TIME_INIT();

    A.ia.reserve(N+1);
    A.ja.reserve(7*N);
    A.a .reserve(7*N);

#define ADD(di,dj,dk,axis) { \
    i1 = index(i+di, j+dj, k+dk); \
    A.ja.push_back(i1); \
    v = 2/(1/k##axis[index_k(i,j,k)] + 1/k##axis[index_k(i+di,j+dj,k+dk)]) / (h##axis * h##axis); \
    A.a.push_back(-v); \
    A.a[dind] += v; \
}

    TIME_START();
    A.nrow = A.ncol = N;
    A.ia.push_back(0);
#ifdef XYZ
    for (uint k = 0; k < nz; k++)
        for (uint j = 0; j < ny; j++)
            for (uint i = 0; i < nx; i++) {
#elif defined ZXY
    for (uint j = 0; j < ny; j++)
        for (uint i = 0; i < nx; i++)
            for (uint k = 0; k < nz; k++) {
#endif
                i0 = index(i, j, k);

                A.ja.push_back(i0);
                A.a .push_back(c);

                uint dind = A.a.size() - 1;
                double v;

#ifdef XYZ
                if (k)	      ADD( 0,  0, -1, z);
                if (j)	      ADD( 0, -1,  0, y);
                if (i)	      ADD(-1,  0,  0, x);
                if (i < nx-1) ADD(+1,  0,  0, x);
                if (j < ny-1) ADD( 0, +1,  0, y);
                if (k < nz-1) ADD( 0,  0, +1, z);
#elif defined ZXY
                if (j)	      ADD( 0, -1,  0, y);
                if (i)	      ADD(-1,  0,  0, x);
                if (k)	      ADD( 0,  0, -1, z);
                if (k < nz-1) ADD( 0,  0, +1, z);
                if (i < nx-1) ADD(+1,  0,  0, x);
                if (j < ny-1) ADD( 0, +1,  0, y);
#endif

                A.ia.push_back(A.ja.size());
            }


    LOG_DEBUG(TIME_INFO("Constructing matrix"));
    LEAVE_MESSAGE("Matrix constructed");
}

/* Construct nonsymmetric M-matrix with diagonal domination
 * Parameter tau's role may vary
 *  - it may serve to determine the random shift percentage in the matrix coefficients
 *  - it may serve as a value q in the speed of convergence (see mini-report on regular splitting)
 */
void SPEMesh::construct_matrix_nonsym(SkylineMatrix& A, double c, double tau) const {
    construct_matrix(A, c);

    LOG_INFO("Constructing nonsymmetric part");

    if (!tau)
        THROW_EXCEPTION("Nonsymmetricity parameter is 0");

    const uvector<uint>& ia = A.ia;
    const uvector<uint>& ja = A.ja;
    uvector<double>&	  a = A.a;

    const uint n = A.size();
    /* Add nonsymmetric part */
    for (uint i = 0; i < n; i++) {
        uint dind = ia[i];

#if 1
        /* Modify each element by small percent */
        for (uint j = ia[i]+1; j < ia[i+1]; j++) {
            double d = random(-tau,tau) * a[j];
            a[j]    -= d;
            a[dind] += d;
        }
#else
        double s = 0.0, c = 0.0;
        /* s = \sum_{j > i} a_{ij} */
        for (uint j = ia[i]; j < ia[i+1]; j++) {
            c += a[j];
            if (ja[j] > i)
                s += (-a[j]);
        }

        if (is_equal(s, 0.0))
            continue;

        double alpha = tau/(1-tau) * c/s;

        a[dind] += alpha * s;
        for (uint j = ia[i]+1; j < ia[i+1]; j++)
            if (ja[j] > i)
                a[j] *= (1 + alpha);
#endif
    }
}

std::ostream& operator<<(std::ostream& os, const Point& p) {
    return os << "(" << p.x << "," << p.y << "," << p.z << ")";
}
