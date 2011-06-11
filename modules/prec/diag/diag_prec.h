#ifndef __DIAG_PREC_H__
#define __DIAG_PREC_H__

#include "modules/prec/prec_base.h"

class DiagPrec: public PrecBase {
private:
    uint n;
    uvector<double> d;	/* d containce inverse of the diagonal */

public:
    DiagPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) const THROW;
};

#endif // __DIAG_PREC_H__
