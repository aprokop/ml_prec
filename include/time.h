#ifndef __TIME_H__
#define __TIME_H__

#include <string>
#include <iomanip>
#include <ctime>

#define TIME_INIT()  clock_t timer
#define TIME_START() timer = clock()
#define TIME_INFO(msg) "TIME: " << msg << " : " << \
	std::fixed << std::setprecision(2) << double(clock() - timer)/CLOCKS_PER_SEC

#endif // __TIME_H__
