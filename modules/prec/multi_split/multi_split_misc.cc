#include "multi_split_prec.h"
#include "modules/prec/misc/misc.h"

#include <map>

void MultiSplitPrec::graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const {
    // construct reverse map
    THROW_EXCEPTION("Need to move to a new graph_planes");
#if 0
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

	::graph_planes(filename, levels[level].A, rev_map, plane, level, mesh);
    } else {
	uint n = level0_A.size();
	for (uint k = 0; k < n; k++)
	    rev_map[k] = k;

	::graph_planes(filename, level0_A, rev_map, plane, level, mesh);
    }
#endif
}

std::ostream& operator<<(std::ostream& os, const MultiSplitPrec& p) {
    os << std::endl;
    os << "nlevels = " << p.nlevels << std::endl;
    float oc = 0, gc = 0;
    for (uint level = 0; level < p.nlevels; level++) {
	gc += p.levels[level].N;
	oc += p.levels[level].nnz;
    }
    gc /= p.levels[0].N;
    oc /= p.levels[0].nnz;
    os << "gc     = " << gc << std::endl;   /* Grid complexity */
    os << "oc     = " << oc << std::endl;   /* Operator complexity */
    for (uint level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	os << p.levels[level];
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const MultiSplitPrec::Level& li) {
    os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
    os << "tr = " << li.tr.size() << ", dtr = " << li.dtr.size() << std::endl;
    os << "q = " << li.q << std::endl;
    os << "ntails = " << li.tails.size() << std::endl;
    os << std::endl;
#if 0
    os << "A: " << li.A;
#endif

    return os;
}
