#include "multi_split_prec.h"

#include <map>

std::ostream& operator<<(std::ostream& os, const MultiSplitPrec& p) {
    os << "nlevels = " << p.nlevels << std::endl;
    for (uint level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	os << p.levels[level];
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const MultiSplitPrec::Level& li) {
    os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
    os << "tr: " << li.tr.size() << ", dtr: " << li.dtr.size() << std::endl;
    os << "q = " << li.q << std::endl;
    os << std::endl;
#if 0
    os << "A: " << li.A;
#endif

    return os;
}
