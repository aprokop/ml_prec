#ifndef __BGS_PREC_H__
#define __BGS_PREC_H__

#include "modules/prec/prec_base.h"

// #define PREC_SUBST

#ifdef PREC_SUBST
#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/prec/diag/diag_prec.h"
#endif

class BGSPrec: public PrecBase {
private:
    uint n, n1, n2;
    SkylineMatrix A1, A2;
    CSRMatrix A21;

#ifdef PREC_SUBST
    PrecBase *B1, *B2;
#endif

    /* UMFPACK factorization of the diagonal blocks */
    mutable void *A1_symbolic, *A1_numeric;
    mutable void *A2_symbolic, *A2_numeric;

public:
    BGSPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) const THROW;
};

#endif // __GS_PREC_H__
