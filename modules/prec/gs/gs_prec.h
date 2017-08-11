#ifndef __GS_PREC_H__
#define __GS_PREC_H__

#include "modules/prec/prec_base.h"

class GSPrec: public PrecBase {
private:
    uint n;
    SkylineMatrix L;

public:
    GSPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) const;
};

#endif // __GS_PREC_H__
