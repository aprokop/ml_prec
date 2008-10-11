#ifndef __MAIN_H__
#define __MAIN_H__

#include "include/define.h"
#include <iostream>
#include <vector>

int set_params(int argc, char * argv[], double& c, double& sigma, 
	       std::vector<double>& sigmas, uint& niter, uint& nwells);

#endif // __MAIN_H__
