#include "include/exception.h"
#include "misc.h"

#include <fstream>
#include <iomanip>

static std::string header(void);
static void mark_set_nodes(const MeshBase& mesh, char plane, uint index, const uvector<uint>& map, uvector<uint>& marked);
void graph_planes(const std::string& filename, const SkylineMatrix& A, const uvector<uint>& map,
		  char plane, const MeshBase& mesh) {
    ASSERT(plane == 'x' || plane == 'y' || plane == 'z', "Unknown plane: " << plane);

    uint n = 0;
    double min_x = 0, max_x = 0, min_y = 0, max_y = 0;;
    switch (plane) {
	case 'x': max_x = mesh.get_size('y'); max_y = mesh.get_size('z'); n = mesh.get_n('x'); break;
	case 'y': max_x = mesh.get_size('x'); max_y = mesh.get_size('z'); n = mesh.get_n('y'); break;
	case 'z': max_x = mesh.get_size('x'); max_y = mesh.get_size('y'); n = mesh.get_n('z'); break;
    }

    double mult_x = 550./(max_x - min_x), mult_y = 690./(max_y - min_y);

#if 0
    /* Preserve scale */
    mult_x = mult_y = std::min(mult_x, mult_y);
#endif

    std::ofstream ofs(filename.c_str());
    if (!ofs.good())
	THROW_EXCEPTION("Cannot open file \"" << filename << "\"");
    ofs << std::fixed << std::setprecision(4);

    ofs << header();

    const std::vector<Point>& nodes = mesh.get_nodes();

    const uvector<uint>& ia = A.get_ia();
    const uvector<uint>& ja = A.get_ja();

    uint pages = 0;

    uint N = A.size();

// #define REGION_ONLY
#ifndef REGION_ONLY
    for (uint k = 0; k < n3; k += 7) {
#else
    uint k = 63; {
#endif

	/* Start page */
	ofs << "%%Page: " << ++pages << std::endl;
	// ofs << "0 setlinewidth\n";
	ofs << "30 40 translate\n";
	ofs << "gsave\n";
	ofs << mult_x << " " << mult_y << " scale\n";

	/* Construct marked array */
	uvector<uint> marked(N, 0);
	mark_set_nodes(mesh, plane, k, map, marked);

	uvector<uint> orth(N, 0);
	for (uint i = 0; i < N; i++)
	    if (marked[i]) {
		const Point& a = nodes[map[i]];

		for (uint j_ = ia[i]+1; j_ < ia[i+1]; j_++) {
		    uint j = ja[j_];
		    if (marked[j]) {
			std::string al = "v";
			if (!A.exist(j,i))
			    al += "v";

			const Point& b = nodes[map[j]];
			switch (plane) {
			    case 'x': ofs << a.y << " " << a.z << " " << b.y << " " << b.z << " " << al << std::endl; break;
			    case 'y': ofs << a.x << " " << a.z << " " << b.x << " " << b.z << " " << al << std::endl; break;
			    case 'z': ofs << a.x << " " << a.y << " " << b.x << " " << b.y << " " << al << std::endl; break;
			}
		    } else {
			orth[i]++;
		    }
		}
	    }

	/* Need to return scale to (1,1) as we'll have ellipses using /a if not */
	ofs << "grestore\n";
#if 0
	for (uint i = 0; i < N; i++)
	    if (marked[i]) {
		switch (orth[i]) {
#if 0
		    case 0: ofs << "b "; break;
#else
		    case 0: continue;
#endif
		    case 1: ofs << "g "; break;
		    case 2: ofs << "r "; break;
		}

		const Point& a = nodes[map[i]];
		switch (plane) {
		    case 'x': ofs << mult_x*a.y << " " << mult_y*a.z << " a\n"; break;
		    case 'y': ofs << mult_x*a.x << " " << mult_y*a.z << " a\n"; break;
		    case 'z': ofs << mult_x*a.x << " " << mult_y*a.y << " a\n"; break;
		}
	    }
#endif
    }

    /* Add footer */
    ofs << "userdict/end-hook known{end-hook}if\n";
    ofs << "%%Pages: " << pages << std::endl;
    ofs << "%%EOF\n";
}

void mark_set_nodes(const MeshBase& mesh, char plane, uint ind, const uvector<uint>& map, uvector<uint>& marked) {
    uint n = map.size();
    uint nx = mesh.get_n('x'), ny = mesh.get_n('y'), nz = mesh.get_n('z');

    std::map<uint,uint> rmap;
    for (uint i = 0; i < n; i++)
	rmap[map[i]] = i;

#ifdef REGION_ONLY
    uint n1min = 44, n1max = 60;
    uint n2min = 50, n2max = 100;
#endif

    std::map<uint,uint>::const_iterator it;
    switch (plane) {
	case 'x':
	    for (uint j = 0; j < ny; j++)
		for (uint k = 0; k < nz; k++)
		    if ((it = rmap.find(mesh.index(ind, j, k))) != rmap.end())
			marked[it->second] = 1;
	    break;
	case 'y':
	    for (uint i = 0; i < nx; i++)
		for (uint k = 0; k < nz; k++)
		    if ((it = rmap.find(mesh.index(i, ind, k))) != rmap.end())
			marked[it->second] = 1;
	    break;
	case 'z':
	    for (uint i = 0; i < nx; i++)
		for (uint j = 0; j < ny; j++) {
#ifdef REGION_ONLY
		    if (i < n1min || i > n1max || j < n2min || j > n2max)
			continue;
#endif
		    if ((it = rmap.find(mesh.index(i, j, ind))) != rmap.end())
			marked[it->second] = 1;
		}
	    break;
    }
}

std::string header(void) {
    std::ostringstream os;

    os << "%!PS-Adobe-2.0\n";
    os << "%%Title: Graph\n";
    os << "%%Creator: prok\n";
    time_t t = time(NULL);
    os << "%%Creation date: " << ctime(&t);
    os << "%%Pages: (atend)\n";
    os << "%%EndComments\n\n";

    os << "%%BeginProlog\n";
    os << "/Helvetica findfont 14 scalefont setfont\n";
    os << "/v {moveto lineto stroke} def\n";
    os << "/ms {moveto show} bind def\n";
    os << "/a{1.1 0 360 arc fill}def\n";
    os << "/r{1 0 0 setrgbcolor}def\n";
    os << "/g{0 1 0 setrgbcolor}def\n";
    os << "/b{0 0 0 setrgbcolor}def\n";
    os << "%%BeginProcSet: arrows 1.0 0\n";
    os << "% \"arrowhead\" takes these arguments:\n";
    os << "% lineweight prevX prevY\n";
    os << "/arrowhead { %def\n";
    os << "gsave\n";
    os << "currentpoint\n";
    os << "4 2 roll exch 4 -1 roll exch\n";
    os << "sub 3 1 roll sub\n";
    os << "exch atan rotate dup scale\n";
    os << "-5 3 rlineto\n";
    os << "1 -3 rlineto\n";
    os << "-1 -3 rlineto\n";
    os << "closepath fill\n";
    os << "grestore\n";
    os << "newpath\n";
    os << "} bind def\n";
    os << "/l^ { %def % lineto-arrow\n";
    os << "currentlinewidth currentpoint 5 3 roll\n";
    os << "lineto\n";
    os << "currentpoint stroke moveto\n";
    os << "arrowhead\n";
    os << "} bind def\n";
    os << "/vv{\n";
    os << "    4 2 roll\n";
    os << "    moveto\n";
    os << "    l^\n";
    os << "} def\n";
    os << "%%EndProcSet: arrows 1.0 0\n";
    os << "%%EndProlog\n\n";
    os << "userdict/start-hook known{start-hook}if\n";

    return os.str();
}
