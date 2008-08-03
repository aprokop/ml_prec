#ifndef __PREC_H__
#define __PREC_H__

#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"
#include <iostream>

class PrecBase {
public:
    virtual ~PrecBase() {};
    virtual void solve(const Vector& f, Vector& x) THROW = 0;
};

class Prec : public PrecBase {
private:
    struct Level {
	uint N;

	Vector x1, u0, u1, f1;
	SkylineMatrix A; // is not set for level 0

	std::vector<int> tr; // local->global
	std::vector<int> dtr;

	double alpha, beta;
	double lmin, lmax;
	uint ncheb;
    };
    uint nlevels;
    std::vector<Level> levels;
    double c;

    double  cheb(double x, uint k) const;
    void    construct_level(uint i, const SkylineMatrix& A);
    void    solve(const Vector& f, Vector& x, uint level) THROW;

public:
    Prec(uint nlevels, double eps, uint ncheb, double c, const SkylineMatrix& A);

    virtual void solve(const Vector& f, Vector& x) THROW;

    friend std::ostream& operator<<(std::ostream& os, const Prec& p);
};

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
    } amg_config;

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
    AMGPrec(const SparseMatrix& A);
    void solve(Vector& f, Vector& x) THROW;
};

#endif
