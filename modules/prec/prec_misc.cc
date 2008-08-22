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

void Prec::graph_planes(uint level, const Mesh& mesh, char plane) const {
    ASSERT(level, "Does not work for initial level");
    ASSERT(plane == 'x' || plane == 'y' || plane == 'z', "Unknown plane: " << plane);

    uint n1 = 0, n2 = 0, n3 = 0;
    double min_x = 0, max_x = 0;
    double min_y = 0, max_y = 0;
    switch (plane) {
	case 'z': n1 = mesh.nx; n2 = mesh.ny; n3 = mesh.nz; 
		  max_x = mesh.size_x; max_y = mesh.size_y;
		  break;
	case 'x': n1 = mesh.ny; n2 = mesh.nz; n3 = mesh.nx;
		  max_x = mesh.size_y; max_y = mesh.size_z;
		  break;
	case 'y': n1 = mesh.nx; n2 = mesh.nz; n3 = mesh.ny; 
		  max_x = mesh.size_x; max_y = mesh.size_z;
		  break;
    }
    double mult_x = 550./(max_x - min_x), mult_y = 690./(max_y - min_y);
#ifdef PRESERVE_SCALE
    mult_x = mult_y = std::min(mult_x, mult_y);
#endif

    std::ofstream ofs("graph.ps");
    ASSERT(ofs.good(), "Cannot open file");
    ofs << std::fixed << std::setprecision(4);

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
    ofs << "/a{1.1 0 360 arc fill}def\n";
    ofs << "/r{1 0 0 setrgbcolor}def\n";
    ofs << "/g{0 1 0 setrgbcolor}def\n";
    ofs << "/b{0 0 0 setrgbcolor}def\n";
    ofs << "%%EndProlog\n\n";

    ofs << "userdict/start-hook known{start-hook}if\n"; 

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

    uint i0 = -1, i1 = -1;
    uint li0, li1;
    std::map<uint,uint>::const_iterator it0, it1;
    const SkylineMatrix& A = levels[level].A;
    ASSERT(A.size() == n, "");
    
    uint left = 0, total = 0;
    uint pages = 0;
    for (uint k = 0; k < n3; k += 7) {
	ofs << "%%Page: " << ++pages << std::endl;
	ofs << "0 setlinewidth\n";
	ofs << "30 40 translate\n";
	ofs << "gsave\n";
	ofs << mult_x << " " << mult_y << " scale\n";

	uint ltotal = 0, lleft = 0;
	for (uint l = 0; l < 2; l++) 
	    for (uint j = 0; j < n2 - l; j++)
		for (uint i = 0; i < n1 - (1-l); i++) {
		    ltotal++;

		    switch (plane) {
			case 'z':
			    i0 = mesh.index(    i,   j, k);
			    i1 = mesh.index(i+1-l, j+l, k);
			    break;
			case 'x':
			    i0 = mesh.index(k,     i,   j);
			    i1 = mesh.index(k, i+1-l, j+l);
			    break;
			case 'y':
			    i0 = mesh.index(    i, k,   j);
			    i1 = mesh.index(i+1-l, k, j+l);
			    break;
		    }

		    if ((it0 = rev_map.find(i0)) == rev_map.end() || 
			(it1 = rev_map.find(i1)) == rev_map.end())
			continue;

		    li0 = it0->second;
		    li1 = it1->second;

		    if (!A.exist(li0, li1))
			continue;

		    lleft++;

		    const Point& a = mesh.nodes[i0], b = mesh.nodes[i1];
		    switch (plane) {
			case 'z': ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " v" << std::endl; break;
			case 'x': ofs << a.y << " " << a.z << " " << b.y << " " << b.z << " v" << std::endl; break;
			case 'y': ofs << a.x << " " << a.z << " " << b.x << " " << b.z << " v" << std::endl; break;
		    }
		}

	// Need to return scale to (1,1) as we'll have ellipses using /a if not
	ofs << "grestore\n";
	for (uint j = 0; j < n2; j++)
	     for (uint i = 0; i < n1; i++) {
		 char fz = 0;
		 switch (plane) {
		     case 'z': i0 = mesh.index(i,j,k); break;
		     case 'x': i0 = mesh.index(k,i,j); break;
		     case 'y': i0 = mesh.index(i,k,j); break;
		 }

		 if ((it0 = rev_map.find(i0)) == rev_map.end())
		     continue;
		 li0 = it0->second;

		 if (k) { 
		     switch (plane) {
			 case 'z': i1 = mesh.index(i,j,k-1); break;
			 case 'x': i1 = mesh.index(k-1,i,j); break;
			 case 'y': i1 = mesh.index(i,k-1,j); break;
		     }
		     if ((it1 = rev_map.find(i1)) != rev_map.end()) {
			 li1 = it1->second;
			 if (A.exist(li0, li1))
			     fz++;
		     }
		 }
		 if (k < n3-1) {
		     switch (plane) {
			 case 'z': i1 = mesh.index(i,j,k+1); break;
			 case 'x': i1 = mesh.index(k+1,i,j); break;
			 case 'y': i1 = mesh.index(i,k+1,j); break;
		     }
		     if ((it1 = rev_map.find(i1)) != rev_map.end()) {
			 li1 = it1->second;
			 if (A.exist(li0, li1))
			     fz++;
		     }
		 }

		 switch (fz) {
		     case 0: continue;
		     case 1: ofs << "g "; break;
		     case 2: ofs << "r "; break;
		 }
		 const Point& a = mesh.nodes[i0];
		 switch (plane) {
		     case 'z': ofs << mult_x*a.x << " " << mult_y*a.y << " a\n"; break;
		     case 'x': ofs << mult_x*a.y << " " << mult_y*a.z << " a\n"; break;
		     case 'y': ofs << mult_x*a.x << " " << mult_y*a.z << " a\n"; break;
		 }
	     }

	ofs << "b (" << plane << " plane #" << k << ", #links = " << lleft << "/" << ltotal << " = " << 
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
