#ifndef __MAIN_H__
#define __MAIN_H__

#include "include/define.h"
#include "modules/common/common.h"
#include <iostream>
#include <vector>

// #define DIAG_PREC
// #define AMG_PREC
#define CHEB_PREC
// #define RELX_PREC

int set_params(int argc, char * argv[], Config& config);

#endif // __MAIN_H__
