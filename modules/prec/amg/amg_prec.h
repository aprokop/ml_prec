#ifndef __AMG_PREC_H__
#define __AMG_PREC_H__

#include "modules/prec/prec_base.h"

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
    AMGConfig amg_config;

    // ***************************************************************
    // Note: we use FORTRAN index enumeration
    // Matrix is stored in skyline format (CSR with diagonal first)
    // ***************************************************************
    std::vector<int> ia, ja;
    std::vector<int> ig;
    std::vector<double> a;

    int n, nnz;

    int nda; 
    int ndia;
    int ndja;
    int ndu;
    int ndf;
    int ndig;

public:
    AMGPrec(const SkylineMatrix& A);

    void solve(Vector& f, Vector& x) THROW;
};

#endif // __AMG_PREC_H__
