#ifndef __TIME_H__
#define __TIME_H__

#include <string>
#include <iomanip>
#include <ctime>

#define TIME_INIT()  clock_t timer
#define TIME_START() timer = clock()
#define TIME_INFO(msg) "TIME: " << msg << " : " << \
        std::fixed << std::setprecision(3) << double(clock() - timer)/CLOCKS_PER_SEC

/* High precision clock (up to nanoseconds) */
inline double pclock() {
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);

    return ts.tv_sec + ts.tv_nsec*1e-9;
}

inline double timediff(const timespec& t1, const timespec& t2) {
    return t2.tv_sec - t1.tv_sec + (t2.tv_nsec - t1.tv_nsec)*1e-9;
}

inline double timediff(const clock_t& t1, const clock_t& t2) {
    return double(t2 - t1)/CLOCKS_PER_SEC;
}

/* More fine-grained macros */
#define TIME_INIT_TIMER(timer) clock_t timer
#define TIME_START_TIMER(timer) timer = clock()
#define TIME_INFO_TIMER(timer,msg) "TIME: " << msg << " : " << \
        std::fixed << std::setprecision(3) << double(clock() - timer)/CLOCKS_PER_SEC

#endif // __TIME_H__
