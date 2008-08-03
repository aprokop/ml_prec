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
		v = 2/(1/kx[i0] + 1/kx[i1]) / (hx*hx);
		A.new_link(i0, i1, v);
	    }
    // y links
    for (int k = 0; k < nz; k++)
	for (int j = 0; j < ny-1; j++)
	    for (int i = 0; i < nx; i++) {
		i0 = index(i, j  , k);
		i1 = index(i, j+1, k);
		v = 2/(1/ky[i0] + 1/ky[i1]) / (hy*hy);
		A.new_link(i0, i1, v);
	    }
    // z links
    for (int k = 0; k < nz-1; k++)
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx; i++) {
		i0 = index(i, j, k);
		i1 = index(i, j, k+1);
		v = 2/(1/kz[i0] + 1/kz[i1]) / (hz*hz);
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
    for (int k = 0; k < nz; k += 7) {
	// ofs << "%%Page: " << k+1 << " " << nz << std::endl;
	ofs << "%%Page: " << ++pages << std::endl;
	ofs << "30 40 translate\n";
	ofs << "gsave\n";
	ofs << mult << " " << mult << " scale\n";

	uint i0, i1;
	uint ltotal = 0, lleft = 0;
	// first direction
	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx-1; i++) {
		ltotal++;

		i0 = index(i, j, k);
		i1 = index(i+1, j, k);

		if (is_removed(i0, i1))
		    continue;
		lleft++;

		const Point& a = nodes[i0], b = nodes[i1];
		ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl;
	    }
	// second direction
	for (int i = 0; i < nx; i++) 
	    for (int j = 0; j < ny-1; j++) {
		ltotal++;

		i0 = index(i, j, k);
		i1 = index(i, j+1, k);

		if (is_removed(i0, i1))
		    continue;
		lleft++;

		const Point& a = nodes[i0], b = nodes[i1];
		ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl;
	    }
	for (int j = 0; j < ny; j++)
	     for (int i = 0; i < nx; i++) {
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
    for (int j = 0; j < ny; j += 10) {
	ofs << "%%Page: " << ++pages << std::endl;
	ofs << "30 40 translate" << std::endl;
	ofs << "gsave\n";
	// we need to scale the problem as size_z << size_x
	ofs << 560./size_x << " " << 720./size_z << " scale\n";
	ofs << "0 setlinewidth\n";

	uint i0, i1;
	uint ltotal = 0, lleft = 0;
	// first direction
	for (int k = 0; k < nz-1; k++)
	    for (int i = 0; i < nx; i++) {
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
