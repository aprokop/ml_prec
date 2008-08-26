#include "config/config.h"
#include "include/logger.h"
#include "include/time.h"
#include "include/tools.h"
#include "mesh.h"

#include <fstream>
#include <algorithm>
#ifdef HAVE_BOOST
#include <boost/lambda/lambda.hpp>
#endif

DEFINE_LOGGER("Mesh");

Mesh::Mesh(double _c) {
    std::ifstream spe("spe_perm.dat");
    ASSERT(spe.good(), "Could not open spe");

    c = _c;

    TIME_INIT();
    TIME_START();
    kx.resize(N);
    ky.resize(N);
    kz.resize(N);
#ifdef HAVE_BOOST
    std::for_each(kx.begin(), kx.end(), spe >> boost::lambda::_1);
    std::for_each(ky.begin(), ky.end(), spe >> boost::lambda::_1);
    std::for_each(kz.begin(), kz.end(), spe >> boost::lambda::_1);
#else
    for (uint i = 0; i < N; i++)
	spe >> kx[i];
    for (uint i = 0; i < N; i++)
	spe >> ky[i];
    for (uint i = 0; i < N; i++)
	spe >> kz[i];
#endif
    LOG_DEBUG(TIME_INFO("Reading SPE file"));
    LEAVE_MESSAGE("Mesh read");

    // nodes
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

void Mesh::construct_matrix(SkylineMatrix& A, uint nwells) const {
    uint i0, i1;
    double v;

    TIME_INIT();

    A.ia.reserve(N+1);
    A.ja.reserve(7*N);
    A.a.reserve(7*N);
    TIME_START();
    A.nrow = A.ncol = N;
    A.ia.push_back(0);
    for (uint j = 0; j < ny; j++)
	for (uint i = 0; i < nx; i++)
	    for (uint k = 0; k < nz; k++) {
		i0 = index(i, j, k);
		A.ja.push_back(i0);
		A.a.push_back(c);

		uint dind = A.a.size() - 1;
		double v;

#define ADD(di,dj,dk,axis) { \
    i1 = index(i+di, j+dj, k+dk); \
    A.ja.push_back(i1); \
    v = 2/(1/k##axis[index_k(i,j,k)] + 1/k##axis[index_k(i+di,j+dj,k+dk)]) / (h##axis * h##axis); \
    A.a.push_back(-v); \
    A.a[dind] += v; \
}
		if (j)	      ADD( 0, -1,  0, y);
		if (i)	      ADD(-1,  0,  0, x);
		if (k)	      ADD( 0,  0, -1, z);
		if (k < nz-1) ADD( 0,  0, +1, z);
		if (i < nx-1) ADD(+1,  0,  0, x);
		if (j < ny-1) ADD( 0, +1,  0, y);

		A.ia.push_back(A.ja.size());
#undef ADD
	    }
    LOG_DEBUG(TIME_INFO("Constructing matrix"));
    LEAVE_MESSAGE("Matrix constructed");

    // Wells
    TIME_START();
    uint nwp = 5;
    v = 100;
    std::vector<uint> ind(nwp);
    for (uint i = 0; i < nwells; i++) {
	LOG_DEBUG("== WELL " << i << " ==");
	uint ix = random(0,nx-1);
	uint iy = random(0,ny-1);
	uint iz = random(0,20);
	ind[0] = index(ix, iy, iz);
	LOG_DEBUG(" (0): " << ix << " " << iy << " " << iz);
	for (uint j = 1; j < nwp; j++) {
	    ix = random(0,nx-1);
	    iy = random(0,ny-1);
	    iz = random(iz, nz-nwp+j);
	    ind[j] = index(ix, iy, iz);
	    LOG_DEBUG(" (" << j << "): " << ix << " " << iy << " " << iz);
	}
	for (uint j = 0; j < nwp-1; j++) {
	    A.add(ind[j], ind[j], v);
	    A.add(ind[j], ind[j+1], -v);
	    A.add(ind[j+1], ind[j], -v);
	    A.add(ind[j+1], ind[j+1], v);
	}
    }
    LOG_DEBUG(TIME_INFO("Adding wells"));
}
