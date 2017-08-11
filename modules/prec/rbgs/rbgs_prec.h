#ifndef __RBGS_PREC_H__
#define __RBGS_PREC_H__

#include "modules/prec/prec_base.h"
#include "project/config.h"

// #define PREC_SUBST

#ifdef PREC_SUBST
#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/prec/diag/diag_prec.h"
#endif

/*
 * Relaxed block Gauss-Seidel preconditioner
 * Same assumptions as for block Gauss-Seidel preconditioner (file bgs_prec.h)
 * The difference is that we allow some elements in the upper triangular block
 */
class RBGSPrec: public PrecBase {
private:
    uint n;
    SkylineMatrix B;

    uint lN;

    static constexpr double alpha = 0.1;

#ifdef PREC_SUBST
    PrecBase *Bprec;
#endif

    /* UMFPACK factorization of the diagonal blocks */
    mutable void *B_symbolic, *B_numeric;

public:
    RBGSPrec(const SkylineMatrix& A, const Config& cfg);

    void solve(Vector& f, Vector& x) const;
};

#endif // __RBGS_PREC_H__
