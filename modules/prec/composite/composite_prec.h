#ifndef __COMPOSITE_PREC_H__
#define __COMPOSITE_PREC_H__

#include "project/config.h"
#include "modules/prec/prec_base.h"

class CompositePrec: public PrecBase {
private:
    uint n;
    std::vector<PrecBase*> precs;

public:
    CompositePrec(const SkylineMatrix& A, const Config& cfg);

    void solve(Vector& f, Vector& x) const;
};

#endif // __COMPOSITE_PREC_H__
