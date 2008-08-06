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

Mesh::Mesh(double _c, uint nwells) {
    std::ifstream spe("spe_perm.dat");
    ASSERT(spe.good(), "Could not open spe");

    c = _c;

    std::vector<double> kx(N), ky(N), kz(N);

    TIME_INIT();

    TIME_START();
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

    uint i0, i1;
    double v;

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
    LEAVE_MESSAGE("Matrix constructed");
    kx.clear();
    ky.clear();
    kz.clear();

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
    LOG_DEBUG(TIME_INFO("Constructing matrix"));
    LEAVE_MESSAGE("Nodes constructed");

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

void Mesh::graph3D() const {
    std::ofstream ofs("graph.gmv");
    ASSERT(ofs.good(), "Cannot open file");

    ofs << "gmvinput ascii" << std::endl << std::endl;

    ofs << "nodev " << nodes.size() << std::endl;
    for (uint i = 0; i < nodes.size(); i++)
	ofs << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z << std::endl;
    ofs << std::endl;

    uint klimit = 2;
    ofs << "faces " << ny*klimit*(nx-1) << " " << ny*klimit*(nx-1) << std::endl;
    // x links
    for (uint k = 0; k < klimit; k++)
	for (uint j = 0; j < ny; j++)
	    for (uint i = 0; i < nx-1; i++) {
		uint i0 = index(i  , j, k) + 1;
		uint i1 = index(i+1, j, k) + 1;
		ofs << "2 " << i0 << " " << i1 << " 1 0" << std::endl;
	    }

    ofs << std::endl << "endgmv" << std::endl;
}

bool Mesh::is_removed(uint i0, uint i1) const {
    ASSERT(i0 != i1, "");

    static double c = 81;
    static double eps = 2;

    if (1 + 12*(-A.get(i0,i1))/c > eps)
	return false;

    return true;
}

void Mesh::graph_xy_planes() const {
    std::ofstream ofs("graph.ps");
    ASSERT(ofs.good(), "Cannot open file");
    ofs << std::fixed << std::setprecision(2);
    ofs << "%!PS-Adobe-2.0\n";
    ofs << "%%Title: Graph\n";
    ofs << "%%Creator: prok\n";
    time_t t = time(NULL);
    ofs << "%%Creation date: " << ctime(&t);
    ofs << "%%Pages: (atend)\n";
    ofs << "%%EndComments\n\n";

    ofs << "%%BeginProlog\n";
    ofs << "/Helvetica findfont 14 scalefont setfont\n";
    ofs << "/v {moveto lineto stroke} def\n";
    ofs << "/ms {moveto show} bind def\n";
    ofs << "/a{3.3 0 360 arc fill}def\n";
    ofs << "/r{1 0 0 setrgbcolor}def\n";
    ofs << "/g{0 1 0 setrgbcolor}def\n";
    ofs << "/b{0 0 0 setrgbcolor}def\n";
    ofs << "%%EndProlog\n\n";

    ofs << "userdict/start-hook known{start-hook}if\n"; 

    // all point are in a rectangle
    double min_x = 0, max_x = size_x;
    double min_y = 0, max_y = size_y;

    ofs << "0 setlinewidth\n";
    double mult = std::min(720./(max_y-min_y), 560./(max_x - min_x));

    
    uint left = 0, total = 0;
    uint pages = 0;
    for (uint k = 0; k < nz; k += 7) {
	// ofs << "%%Page: " << k+1 << " " << nz << std::endl;
	ofs << "%%Page: " << ++pages << std::endl;
	ofs << "30 40 translate\n";
	ofs << "gsave\n";
	ofs << mult << " " << mult << " scale\n";

	uint i0, i1;
	uint ltotal = 0, lleft = 0;
	// first direction
	for (uint j = 0; j < ny; j++)
	    for (uint i = 0; i < nx-1; i++) {
		ltotal++;


		i0 = index(i, j, k);
		i1 = index(i+1, j, k);

		if (k == 0)
		    LOG_DEBUG("(" << i << "," << j << ") - (" << i+1 << "," << j <<") : " << A.get(i0,i1));

		if (is_removed(i0, i1))
		    continue;
		lleft++;

		const Point& a = nodes[i0], b = nodes[i1];
		ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl;
	    }
	// second direction
	for (uint i = 0; i < nx; i++) 
	    for (uint j = 0; j < ny-1; j++) {
		ltotal++;

		i0 = index(i, j, k);
		i1 = index(i, j+1, k);

		if (is_removed(i0, i1))
		    continue;
		lleft++;

		const Point& a = nodes[i0], b = nodes[i1];
		ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl;
	    }
	for (uint j = 0; j < ny; j++)
	     for (uint i = 0; i < nx; i++) {
		 char fz = 0;
		 i0 = index(i,j,k);
		 if (k) {
		     i1 = index(i,j,k-1);
		     if (!is_removed(i0, i1))
			 fz++;
		 }
		 if (k < nz-1) {
		     i1 = index(i,j,k+1);
		     if (!is_removed(i0, i1))
			 fz++;
		 }
		 switch (fz) {
		     case 0: continue;
		     case 1: ofs << "g "; break;
		     case 2: ofs << "r "; break;
		 }
		 const Point& a = nodes[i0];
		 ofs << a.x << " " << a.y << " a\n";
	     }

	ofs << "grestore\n";
	ofs << "b (Level = " << k << ", #links = " << lleft << "/" << ltotal << " = " << 
		100.*lleft/ltotal << "%) 10 730 ms" << std::endl;
	ofs << "showpage" << std::endl;
	ofs << "%%EndPage\n";
	
	total += ltotal;
	left  += lleft;
    }
    ofs << "userdict/end-hook known{end-hook}if\n";
    ofs << "%%Pages: " << pages << std::endl;
    ofs << "%%EOF\n";

    std::cout << "#links = " << left << "/" << total << " = " << 100.*left/total << "%" << std::endl;
}

void Mesh::graph_z_lines() const {
    std::ofstream ofs("graph.ps");
    ASSERT(ofs.good(), "Cannot open file");
    ofs << std::fixed << std::setprecision(2);
    ofs << "%!PS-Adobe-2.0\n";
    ofs << "%%Title: Graph\n";
    ofs << "%%Creator: prok\n";
    time_t t = time(NULL);
    ofs << "%%Creation date: " << ctime(&t);
    ofs << "%%Pages: (atend)\n";
    ofs << "%%EndComments\n\n";

    ofs << "%%BeginProlog\n";
    ofs << "/Helvetica findfont 14 scalefont setfont\n";
    ofs << "/v {moveto lineto stroke} def\n";
    ofs << "/ms {moveto show} bind def\n";
    ofs << "%%EndProlog\n\n";

    ofs << "userdict/start-hook known{start-hook}if\n"; 

    
    uint left = 0, total = 0;
    uint pages = 0;
    for (uint j = 0; j < ny; j += 10) {
	ofs << "%%Page: " << ++pages << std::endl;
	ofs << "30 40 translate" << std::endl;
	ofs << "gsave\n";
	// we need to scale the problem as size_z << size_x
	ofs << 560./size_x << " " << 720./size_z << " scale\n";
	ofs << "0 setlinewidth\n";

	uint i0, i1;
	uint ltotal = 0, lleft = 0;
	// first direction
	for (uint k = 0; k < nz-1; k++)
	    for (uint i = 0; i < nx; i++) {
		ltotal++;

		i0 = index(i, j, k);
		i1 = index(i, j, k+1);

		if (is_removed(i0, i1))
		    continue;
		lleft++;

		const Point& a = nodes[i0], b = nodes[i1];

		ofs << a.x << " " << a.z << " " << b.x << " " << b.z << " v" << std::endl;
	    }

	ofs << "grestore\n";
	ofs << "(j = " << j << ", #links = " << lleft << "/" << ltotal << " = " << 
		100.*lleft/ltotal << "%) 10 720 ms" << std::endl;
	ofs << "showpage" << std::endl;
	ofs << "%%EndPage\n";
	
	total += ltotal;
	left  += lleft;
    }
    ofs << "userdict/end-hook known{end-hook}if\n";
    ofs << "%%Pages: " << pages << std::endl;
    ofs << "%%EOF\n";

    std::cout << "#links = " << left << "/" << total << " = " << 100.*left/total << "%" << std::endl;
}
