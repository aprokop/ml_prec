#include "prec.h"
#include "include/logger.h"

#include <cmath>

DEFINE_LOGGER("Prec");

std::ostream& operator<<(std::ostream& os, const Prec& p) {
    os << "nlevels = " << p.nlevels;
    for (uint level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	const Prec::Level& li = p.levels[level];
	os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
	os << "tr: " << li.tr.size() << ", dtr: " << li.dtr.size() << std::endl;
	if (level != p.nlevels-1) {
	    os << "Ncheb = " << li.ncheb << ", ";
	}
	os << "[lmin, lmax] = [" << li.lmin << "," << li.lmax << "], cond = " << li.lmax/li.lmin << std::endl;
	os << "alpha = " << li.alpha << ", beta = " << li.beta << std::endl;
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
    if (p.levels[p.nlevels-1].ncheb)
	std::cout << "Coarse ncheb = " << p.levels[p.nlevels-1].ncheb << std::endl;
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
