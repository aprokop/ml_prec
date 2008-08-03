#ifndef __PREC_H__
#define __PREC_H__

#include "modules/vector/vector.h"
#include "modules/matrix/matrix.h"
#include <iostream>

class PrecBase {
public:
    virtual ~PrecBase() {};
    virtual void solve(const Vector& f, Vector& x) const THROW = 0;
};

class Prec : public PrecBase {
private:
    struct Level {
	uint N;

	std::vector<double> x1, f1;
	SkylineMatrix& A; // is not set for level 0

	std::vector<int> tr; // local->global
	std::vector<int> dtr;

	double alpha, beta;
	double lmin, lmax;
    };

    void solve(const Vector& f, Vector& x, uint level) const THROW;

public:
    Prec(uint nlevels, double eps, const SkylineMatrix& A);

    virtual void solve(const Vector& f, Vector& x) const THROW;

    friend std::ostream& operator<<(std::ostream& os, const Prec& p);
};

#endif
