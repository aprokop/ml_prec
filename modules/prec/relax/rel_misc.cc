#include "rel_prec.h"
#include "modules/prec/misc/misc.h"
#include "include/logger.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>

DEFINE_LOGGER("Prec");

std::ostream& operator<<(std::ostream& os, const RelPrec& p) {
    os << std::endl;
    os << "nlevels = " << p.nlevels << std::endl;
    os << "niter   = " << p.niter << std::endl;
    os << "gamma   = " << p.gamma << std::endl;
    os << "sigmas  = ";
    for (uint i = 0; i < p.sigmas.size(); i++)
        os << p.sigmas[i] << " ";
    os << std::endl;
    for (uint level = 0; level < p.nlevels; level++) {
        os << std::endl << "================== Level: " << level << " =======================" << std::endl;
        const RelPrec::Level& li = p.levels[level];
        os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
        os << "tr: " << li.tr.size() << ", dtr: " << li.dtr.size() << std::endl;
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

void RelPrec::graph_planes(const std::string& filename, uint level, char plane) const {
    // construct reverse map
    std::map<uint,uint> rev_map;
    if (level) {
        std::vector<uint> gtr = levels[level-1].tr;
        uint n = gtr.size();
        for (int l = level-2; l >= 0; l--)
            for (uint i = 0; i < n; i++)
                gtr[i] = levels[l].tr[gtr[i]];
        for (uint i = 0; i < n; i++)
            rev_map[gtr[i]] = i;
        gtr.clear();
    }

    ::graph_planes(filename, levels[level].A, rev_map, plane, level, mesh);
}
