#ifndef __PREC_BASE_H__
#define __PREC_BASE_H__

#include "include/define.h"
#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"

class PrecBase {
public:
    virtual ~PrecBase() {};
    virtual void solve(Vector& f, Vector& x) const THROW = 0;
};

#endif // __PREC_BASE_H__
