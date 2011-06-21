#ifndef __RBGS_PREC_H__
#define __RBGS_PREC_H__

#include "modules/prec/prec_base.h"

// #define PREC_SUBST

#ifdef PREC_SUBST
#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/prec/diag/diag_prec.h"
#endif

class RBGSPrec: public PrecBase {
private:
    uint n;
    SkylineMatrix B;

    static const double alpha = 0.95;

#ifdef PREC_SUBST
    PrecBase *Bprec;
#endif

    /* UMFPACK factorization of the diagonal blocks */
    mutable void *B_symbolic, *B_numeric;

public:
    RBGSPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) const THROW;
};

#endif // __RBGS_PREC_H__
