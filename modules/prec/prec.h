#ifndef __PREC_H__
#define __PREC_H__

#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"
#include "modules/mesh/mesh.h"
#include <iostream>

class PrecBase {
public:
    virtual ~PrecBase() {};
    virtual void solve(Vector& f, Vector& x) THROW = 0;
};

class Prec : public PrecBase {
private:
    double c;
    double galpha, gbeta;

    struct Level {
	uint N, nnz;

	Vector x1, u0, u1, f1;
	SkylineMatrix A; // is not set for level 0

	std::vector<uint> tr; 
	std::vector<uint> dtr;

	double lmin, lmax;
    };
    uint nlevels;
    std::vector<Level> levels;
    uint ncheb;

    double  cheb(double x, uint k) const;
    void    construct_level(uint i, const SkylineMatrix& A);
    void    solve(const Vector& f, Vector& x, uint level) THROW;

public:
    Prec(double eps, uint ncheb, double c, const SkylineMatrix& A);

    void graph_xy_planes(uint level, const Mesh& mesh, const SkylineMatrix& A) const;
    // void graph_z_lines() const;

    virtual void solve(Vector& f, Vector& x) THROW;

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

    virtual void solve(Vector& f, Vector& x) THROW;
};

#endif
