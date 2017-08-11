#include "include/logger.h"
#include "include/tools.h"
#include "composite_prec.h"

#include "modules/prec/amg/amg_prec.h"
#include "modules/prec/bgs/bgs_prec.h"
#include "modules/prec/cheb/cheb_prec.h"
#include "modules/prec/diag/diag_prec.h"
#include "modules/prec/gs/gs_prec.h"
#include "modules/prec/id/id_prec.h"
#include "modules/prec/multi_split/multi_split_prec.h"
#include "modules/prec/rbgs/rbgs_prec.h"
#include "modules/prec/relax/rel_prec.h"
#include "modules/prec/sym/sym_prec.h"

DEFINE_LOGGER("CompositePrec");

void CompositePrec::solve(Vector& f, Vector& x) const{
    ASSERT(f.size() == n && x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    uvector<double> x1(x.size());
    precs[0]->solve(f, x);
    for (uint i = 1; i < precs.size(); i++) {
        x1.swap(x);
        precs[i]->solve(x1, x);
    }
}

CompositePrec::CompositePrec(const SkylineMatrix& A, const Config& cfg) {
    ASSERT(A.cols() == A.rows(), "Matrix must be square");
    n = A.rows();

    precs.resize(2);
    precs[0] = new GSPrec(A);
    precs[1] = new MultiSplitPrec(A, cfg);
}
