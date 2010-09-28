#include "multi_split_prec.h"
#include "modules/prec/misc/misc.h"

#include <map>

void MultiSplitPrec::graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const {
    if (level >= nlevels)
	THROW_EXCEPTION("Expected level < " << nlevels);

    uvector<uint> map;
    if (level) {
	uint n = levels[level].A.size();
	const Level& lp = levels[level-1];

	/* Construct map: indices on the level -> indices in the original (level 0) matrix */
	map.resize(n);
	memcpy(&map[0], &lp.map[lp.M], (lp.N - lp.M - lp.Md)*sizeof(uint));

	for (int l = level-2; l >= 0; l--) {
	    const Level& li = levels[l];
	    for (uint i = 0; i < n; i++)
		map[i] = li.map[li.M + map[i]];
	}

	/* Call the plotting procedure */
	::graph_planes(filename, levels[level].A, map, plane, mesh);
    } else {
	uint n = level0_A.size();

	/* Construct map */
	map.resize(n);
	for (uint i = 0; i < n; i++)
	    map[i] = i;

	::graph_planes(filename, level0_A, map, plane, mesh);
    }
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
    os << "M = " << li.M << ", Md  = " << li.Md << std::endl;
    os << "q = " << li.q << std::endl;
    os << std::endl;
#if 0
    os << "A: " << li.A;
#endif

    return os;
}
