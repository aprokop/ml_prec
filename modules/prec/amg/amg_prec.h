#ifndef __AMG_PREC_H__
#define __AMG_PREC_H__

#include "modules/prec/prec_base.h"
#include "include/uvector.h"

class AMGPrec: public PrecBase {
private:
    struct AMGConfig {
	int matrix;
	int iswtch;
	int iout;
	int iprint;
	int levelx;
	int ifirst;
	int ncyc;
	double eps;
	int madapt;
	int nrd;
	int nsolco;
	int nru;
	double ecg1;
	double ecg2;
	double ewt2;
	int nwt;
	int ntr;

	AMGConfig();
    };
    mutable AMGConfig amg_config;

    // ***************************************************************
    // Note: we use FORTRAN index enumeration
    // Matrix is stored in skyline format (CSR with diagonal first)
    // ***************************************************************
    mutable uvector<int> ia, ja;
    mutable uvector<int> ig;
    mutable uvector<double> a;

    mutable int n, nnz;

    mutable int nda;
    mutable int ndia;
    mutable int ndja;
    mutable int ndu;
    mutable int ndf;
    mutable int ndig;

public:
    AMGPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) const THROW;
};

#endif // __AMG_PREC_H__
