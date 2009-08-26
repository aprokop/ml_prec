#ifndef __MAIN_H__
#define __MAIN_H__

#include "include/define.h"
#include "modules/common/common.h"
#include <iostream>
#include <vector>

int set_params(int argc, char * argv[], Config& config);
double avg_time(const std::vector<double>& times);

#endif // __MAIN_H__
