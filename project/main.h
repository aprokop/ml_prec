#ifndef __MAIN_H__
#define __MAIN_H__

#include "include/define.h"
#include <iostream>
#include <vector>

// #define DIAG_PREC
#define AMG_PREC
// #define CHEB_PREC
// #define RELX_PREC

struct Config {
    uint   niter; 
    uint   nthreads;
    uint   nwells;
    uint   nx, ny, nz;
    double c; 
    double sigma; 
    double gamma;
    std::vector<double> sigmas;

    friend std::ostream& operator<<(std::ostream& os, const Config& config);
};

int set_params(int argc, char * argv[], Config& config);

#endif // __MAIN_H__
