#include "prec.h"
#include "include/logger.h"

#include <cmath>
#include <fstream>
#include <iomanip>

DEFINE_LOGGER("Prec");

std::ostream& operator<<(std::ostream& os, const Prec& p) {
    os << "nlevels = " << p.nlevels << std::endl;
    os << "alpha = " << p.galpha << ", beta = " << p.gbeta << std::endl;
    os << "Ncheb = " << p.ncheb;
    for (uint level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	const Prec::Level& li = p.levels[level];
	os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
	os << "tr: " << li.tr.size() << ", dtr: " << li.dtr.size() << std::endl;
	os << "[lmin, lmax] = [" << li.lmin << "," << li.lmax << "], cond = " << li.lmax/li.lmin << std::endl;
#if 0
	if (level < p.nlevels-1) {
	    os << "tr: " << li.tr;
	    os << "dtr: " << li.dtr;
	}
#endif
#if 0
	if (level) 
	    os << "A: " << li.A;
#endif
    }
    return os;
} 

double Prec::cheb(double x, uint k) const {
   ASSERT(x >= 1, "");
   // return cosh(k*acosh(x));

   switch(k) {
       case 0:	return 1;
       case 1:	return x;
       case 2:	return 2*x*x-1;
       case 3:	return x*(4*x*x - 3);
       default: return cosh(k*acosh(x));
   }
}

void Prec::graph_xy_planes(uint level, const Mesh& mesh, const SkylineMatrix& B) const {
    ASSERT(level, "Does not work for initial level");

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
    double min_x = 0, max_x = mesh.size_x;
    double min_y = 0, max_y = mesh.size_y;

    ofs << "0 setlinewidth\n";
    double mult = std::min(720./(max_y - min_y), 560./(max_x - min_x));

    // construct reverse map
    std::vector<uint> gtr = levels[level-1].tr;
    uint n = gtr.size();
    for (int l = level-2; l >= 0; l--)
	for (uint i = 0; i < n; i++)
	    gtr[i] = levels[l].tr[gtr[i]];
    std::map<uint,uint> rev_map;
    for (uint i = 0; i < n; i++)
	rev_map[gtr[i]] = i;
    gtr.clear();

    uint i0, i1;
    uint li0, li1;
    std::map<uint,uint>::const_iterator it0, it1;
    const SkylineMatrix& A = levels[level].A;
    ASSERT(A.size() == n, "");
    
    uint left = 0, total = 0;
    uint pages = 0;
    for (uint k = 0; k < mesh.nz; k += 700) {
	ofs << "%%Page: " << ++pages << std::endl;
	ofs << "30 40 translate\n";
	ofs << "gsave\n";
	ofs << mult << " " << mult << " scale\n";

	uint ltotal = 0, lleft = 0;
	// first direction
	for (uint j = 0; j < mesh.ny; j++)
	    for (uint i = 0; i < mesh.nx-1; i++) {
		ltotal++;

		i0 = mesh.index(i,   j, k);
		i1 = mesh.index(i+1, j, k);

		if ((it0 = rev_map.find(i0)) == rev_map.end() || 
		    (it1 = rev_map.find(i1)) == rev_map.end())
		    continue;

		li0 = it0->second;
		li1 = it1->second;

		if (!A.exist(li0, li1))
		    continue;

		lleft++;

		const Point& a = mesh.nodes[i0], b = mesh.nodes[i1];
		ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl;
	    }

	// second direction
	for (uint i = 0; i < mesh.nx; i++) 
	    for (uint j = 0; j < mesh.ny-1; j++) {
		ltotal++;

		i0 = mesh.index(i,   j, k);
		i1 = mesh.index(i, j+1, k);

		if ((it0 = rev_map.find(i0)) == rev_map.end() || 
		    (it1 = rev_map.find(i1)) == rev_map.end())
		    continue;

		li0 = it0->second;
		li1 = it1->second;

		if (!A.exist(li0, li1))
		    continue;

		lleft++;

		const Point& a = mesh.nodes[i0], b = mesh.nodes[i1];
		ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl;
	    }

	for (uint j = 0; j < mesh.ny; j++)
	     for (uint i = 0; i < mesh.nx; i++) {
		 char fz = 0;
		 i0 = mesh.index(i,j,k);
		 if ((it0 = rev_map.find(i0)) == rev_map.end())
		     continue;
		 li0 = it0->second;

		 if (k && (it1 = rev_map.find(mesh.index(i,j,k-1))) != rev_map.end()) {
		     li1 = it1->second;
		     if (A.exist(li0, li1))
			 fz++;
		 }
		 if (k < mesh.nz-1 && (it1 = rev_map.find(mesh.index(i,j,k+1))) != rev_map.end()) {
		     li1 = it1->second;
		     if (A.exist(li0, li1))
			 fz++;
		 }
		 switch (fz) {
		     case 0: continue;
		     case 1: ofs << "g "; break;
		     case 2: ofs << "r "; break;
		 }
		 const Point& a = mesh.nodes[i0];
		 ofs << a.x << " " << a.y << " a\n";
	     }

	ofs << "grestore\n";
	ofs << "b (z level = " << k << ", #links = " << lleft << "/" << ltotal << " = " << 
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

// void Mesh::graph_z_lines() const {
    // std::ofstream ofs("graph.ps");
    // ASSERT(ofs.good(), "Cannot open file");
    // ofs << std::fixed << std::setprecision(2);
    // ofs << "%!PS-Adobe-2.0\n";
    // ofs << "%%Title: Graph\n";
    // ofs << "%%Creator: prok\n";
    // time_t t = time(NULL);
    // ofs << "%%Creation date: " << ctime(&t);
    // ofs << "%%Pages: (atend)\n";
    // ofs << "%%EndComments\n\n";
// 
    // ofs << "%%BeginProlog\n";
    // ofs << "/Helvetica findfont 14 scalefont setfont\n";
    // ofs << "/v {moveto lineto stroke} def\n";
    // ofs << "/ms {moveto show} bind def\n";
    // ofs << "%%EndProlog\n\n";
// 
    // ofs << "userdict/start-hook known{start-hook}if\n"; 
// 
    // 
    // uint left = 0, total = 0;
    // uint pages = 0;
    // for (uint j = 0; j < ny; j += 10) {
	// ofs << "%%Page: " << ++pages << std::endl;
	// ofs << "30 40 translate" << std::endl;
	// ofs << "gsave\n";
	// // we need to scale the problem as size_z << size_x
	// ofs << 560./size_x << " " << 720./size_z << " scale\n";
	// ofs << "0 setlinewidth\n";
// 
	// uint i0, i1;
	// uint ltotal = 0, lleft = 0;
	// // first direction
	// for (uint k = 0; k < nz-1; k++)
	    // for (uint i = 0; i < nx; i++) {
		// ltotal++;
// 
		// i0 = index(i, j, k);
		// i1 = index(i, j, k+1);
// 
		// if (is_removed(i0, i1))
		    // continue;
		// lleft++;
// 
		// const Point& a = nodes[i0], b = nodes[i1];
// 
		// ofs << a.x << " " << a.z << " " << b.x << " " << b.z << " v" << std::endl;
	    // }
// 
	// ofs << "grestore\n";
	// ofs << "(j = " << j << ", #links = " << lleft << "/" << ltotal << " = " << 
		// 100.*lleft/ltotal << "%) 10 720 ms" << std::endl;
	// ofs << "showpage" << std::endl;
	// ofs << "%%EndPage\n";
	// 
	// total += ltotal;
	// left  += lleft;
    // }
    // ofs << "userdict/end-hook known{end-hook}if\n";
    // ofs << "%%Pages: " << pages << std::endl;
    // ofs << "%%EOF\n";
// 
    // std::cout << "#links = " << left << "/" << total << " = " << 100.*left/total << "%" << std::endl;
// }
