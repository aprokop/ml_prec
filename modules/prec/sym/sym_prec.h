#ifndef __SYM_H__
#define __SYM_H__

#include "modules/prec/prec_base.h"
#include "modules/prec/cheb/cheb_prec.h"
#include <memory>

/*
 * Given an unsymmetric M-matrix constructs an easy way to invert its symmetric part
 * with the accuracy eps
 */
class SymPrec: public PrecBase {
private:
    uint		n;
    SkylineMatrix	Asym;
    std::auto_ptr<Prec>	B;

    double		eps;

public:
    SymPrec(const SkylineMatrix& A, const Config& cfg);

    void set_eps(double eps_) {
        eps = eps_;
    }
    void solve(Vector& f, Vector& x) const THROW;
};

#endif // __SYM_H__
