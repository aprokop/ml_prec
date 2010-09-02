#ifndef __DIAG_PREC_H__
#define __DIAG_PREC_H__

#include "modules/prec/prec_base.h"

class DiagPrec: public PrecBase {
private:
    uint n;
    std::vector<double> d;

public:
    DiagPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) const THROW;
};

#endif // __DIAG_PREC_H__
