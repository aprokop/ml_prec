#include "include/logger.h"
#include "include/time.h"
#include "include/tools.h"
#include "mesh.h"

#include <fstream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

DEFINE_LOGGER("Mesh");

Mesh::Mesh() {
    std::ifstream spe("spe_perm.dat");
    ASSERT(spe.good(), "Could not open spe");

    std::vector<double> kx(N), ky(N), kz(N);

    TIME_INIT();

    TIME_START();
    std::for_each(kx.begin(), kx.end(), spe >> boost::lambda::_1);
    std::for_each(ky.begin(), ky.end(), spe >> boost::lambda::_1);
    std::for_each(kz.begin(), kz.end(), spe >> boost::lambda::_1);
    LOG_DEBUG(TIME_INFO("Reading SPE file"));
    LEAVE_MESSAGE("Mesh read");

    uint i0, i1;
    double v;

    TIME_START();
    A = FVSparseMatrix(N);
    LEAVE_MESSAGE("Matrix allocated");
    // x links
    for (int k = 0; k < nz; k++)
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx-1; i++) {
		i0 = index(i  , j, k);
		i1 = index(i+1, j, k);
		v = hy*hz * 2/hx * 1/(1/kx[i0] + 1/kx[i1]);
		A.new_link(i0, i1, v);
	    }
    // y links
    for (int k = 0; k < nz; k++)
	for (int j = 0; j < ny-1; j++)
	    for (int i = 0; i < nx; i++) {
		i0 = index(i, j  , k);
		i1 = index(i, j+1, k);
		v = hx*hz * 2/hy * 1/(1/ky[i0] + 1/ky[i1]);
		A.new_link(i0, i1, v);
	    }
    // z links
    for (int k = 0; k < nz-1; k++)
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx; i++) {
		i0 = index(i, j, k);
		i1 = index(i, j, k+1);
		v = hx*hy * 2/hz * 1/(1/kz[i0] + 1/kz[i1]);
		A.new_link(i0, i1, v);
	    }
    LEAVE_MESSAGE("Matrix constructed");
    // nodes
    nodes.resize(N);
    for (int k = 0; k < nz; k++)
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx; i++) {
		uint ind = index(i,j,k);
		nodes[ind].x = i*hx;
		nodes[ind].y = j*hy;
		nodes[ind].z = k*hz;
	    }
    LOG_DEBUG(TIME_INFO("Constructing matrix"));
    LEAVE_MESSAGE("Nodes constructed");
}

void Mesh::graph3D() const {
    std::ofstream ofs("graph.gmv");
    ASSERT(ofs.good(), "Cannot open file");

    ofs << "gmvinput ascii" << std::endl << std::endl;

    ofs << "nodev " << nodes.size() << std::endl;
    for (int i = 0; i < nodes.size(); i++)
	ofs << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z << std::endl;
    ofs << std::endl;

    int klimit = 2;
    ofs << "faces " << ny*klimit*(nx-1) << " " << ny*klimit*(nx-1) << std::endl;
    // x links
    for (int k = 0; k < klimit; k++)
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx-1; i++) {
		uint i0 = index(i  , j, k) + 1;
		uint i1 = index(i+1, j, k) + 1;
		ofs << "2 " << i0 << " " << i1 << " 1 0" << std::endl;
	    }

    ofs << std::endl << "endgmv" << std::endl;
}

void Mesh::graph_2D_plane(char type, uint pos) const {
    std::ofstream ofs("graph.ps");
    ASSERT(ofs.good(), "Cannot open file");
    ofs << std::fixed << std::setprecision(2);
    ofs << "%!PS-Adobe-2.0\n";
    ofs << "%%BeginProlog\n";
    ofs << "/Helvetica findfont 14 scalefont setfont\n\n";
    ofs << "/v {moveto lineto stroke} def\n";
    ofs << "/ms {moveto show} bind def\n";
    ofs << "%%EndProlog\n";

    ofs << "userdict/start-hook known{start-hook}if\n"; 

    // all point are in a rectangle
    double min_x = 0, max_x = size_x;
    double min_y = 0, max_y = size_y;

    double mult = std::min(800./(max_y-min_y), 560./(max_x - min_x));
    // LOG_DEBUG("mult = " << mult);

    ofs << "20 20 translate" << std::endl;
    ofs << "gsave\n";
    
    for (int k = 0; k < nz; k++) {
	ofs << "%%Page: " << k+1 << " " << nz << std::endl;
	ofs << "gsave\n";
	ofs << "1.0 1.0 scale 0 setlinewidth\n";
	ofs << "0 setlinewidth" << std::endl;

	uint i0, i1;
	uint total = 0, left = 0;
	// first direction
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx-1; i++) {
		total++;

		i0 = index(i, j, k);
		i1 = index(i+1, j, k);

		if (1 - 12*A.get(i0,i1) > 3)
		    continue;
		left++;

		const Point& a = nodes[i0], b = nodes[i1];

		ofs <<  mult * (a.x - min_x) << " " << mult * (a.y - min_y) << " " <<
			mult * (b.x - min_x) << " " << mult * (b.y - min_y) << " v" << std::endl;
	    }
	// second direction
	for (int i = 0; i < nx; i++) 
	    for (int j = 0; j < ny-1; j++) {
		total++;

		i0 = index(i, j, k);
		i1 = index(i, j+1, k);

		if (1 - 12*A.get(i0,i1) > 3)
		    continue;
		left++;

		const Point& a = nodes[i0], b = nodes[i1];

		ofs <<	mult * (a.x - min_x) << " " << mult * (a.y - min_y) << " " <<
			mult * (b.x - min_x) << " " << mult * (b.y - min_y) << " v" << std::endl;
	    }
	ofs << "(Level = " << k << ", #links = " << left << "/" << total << " = " << 
		100.*left/total << "%) 10 810 ms" << std::endl;
	ofs << "showpage" << std::endl;
	ofs << "grestore\n";
	ofs << "%%EndPage\n";
    }
    ofs << "userdict/end-hook known{end-hook}if\n";
    ofs << "%%EOF\n";
}
