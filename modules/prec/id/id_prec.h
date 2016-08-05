#ifndef __ID_PREC_H__
#define __ID_PREC_H__

#include "modules/prec/prec_base.h"

class IdPrec: public PrecBase {
private:
    uint n;

public:
    IdPrec(const SkylineMatrix& A) {
        n = A.size();
    }

    void solve(Vector& f, Vector& x) const THROW {
        ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());
        x = f;
    }
};

#endif // __ID_PREC_H__
