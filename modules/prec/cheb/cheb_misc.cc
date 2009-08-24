#include "cheb_prec.h"
#include "modules/prec/misc/misc.h"
#include "include/logger.h"

#include <map>

DEFINE_LOGGER("Prec");

void Prec::graph_planes(const std::string& filename, uint level, char plane) const {
    // construct reverse map
    std::map<uint,uint> rev_map;
    if (level) {
	uvector<uint> gtr = levels[level-1].tr;
	uint n = gtr.size();
	for (int l = level-2; l >= 0; l--)
	    for (uint i = 0; i < n; i++)
		gtr[i] = levels[l].tr[gtr[i]];
	for (uint i = 0; i < n; i++)
	    rev_map[gtr[i]] = i;
	gtr.clear();
    } 

    // ::graph_planes(filename, levels[level].A, rev_map, plane, level, mesh);
}

std::ostream& operator<<(std::ostream& os, const Prec& p) {
    os << "nlevels = " << p.nlevels << std::endl;
    for (uint level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	os << p.levels[level];
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const Prec::Level& li) {
    os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
    os << "tr: " << li.tr.size() << ", dtr: " << li.dtr.size() << std::endl;
    os << "ncheb = " << li.ncheb << std::endl;
    os << "alpha = " << li.alpha << ", beta = " << li.beta << std::endl;
    os << "[lmin, lmax] = [" << li.lmin << "," << li.lmax << "], cond = " << li.lmax/li.lmin << std::endl;
    os << std::endl;
#if 0
    os << "A: " << li.A;
#endif
    return os;
}
