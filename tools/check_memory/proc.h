#ifndef PROC_H
#define PROC_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/stat.h>

typedef struct proc_t {
    unsigned long long
	utime,		// stat            user-mode CPU time accumulated by process
	stime,		// stat            kernel-mode CPU time accumulated by process
	cutime,		// stat            cumulative utime of process and reaped children
	cstime,		// stat            cumulative stime of process and reaped children
	start_time;	// stat            start time of process -- seconds since 1-1-70
    long
	rss;		// stat            resident set size from /proc/#/stat (pages)
} proc_t;

proc_t * get_proc_stats(pid_t pid, proc_t *p);

#endif
