/* ALL PROCEDURES ARE FROM procps 3.2.7 */

#include "proc.h"
#include <string.h>

typedef struct proc_stat_t {
    int
        tid,		// (special)       task id, the POSIX thread ID (see also: tgid)
    	ppid;		// stat,status     pid of parent process
    unsigned
        pcpu;           // stat (special)  %CPU usage (is not filled in by readproc!!!)
    char
    	state;		// stat,status     single-char code for process state (S=sleeping)
    unsigned long long
	utime,		// stat            user-mode CPU time accumulated by process
	stime,		// stat            kernel-mode CPU time accumulated by process
// and so on...
	cutime,		// stat            cumulative utime of process and reaped children
	cstime,		// stat            cumulative stime of process and reaped children
	start_time;	// stat            start time of process -- seconds since 1-1-70
    unsigned long
	start_code,	// stat            address of beginning of code segment
	end_code,	// stat            address of end of code segment
	start_stack,	// stat            address of the bottom of stack for the process
	kstk_esp,	// stat            kernel stack pointer
	kstk_eip,	// stat            kernel instruction pointer
	wchan;		// stat (special)  address of kernel wait channel proc is sleeping in
    long
	priority,	// stat            kernel scheduling priority
	nice,		// stat            standard unix nice level of process
	rss,		// stat            resident set size from /proc/#/stat (pages)
	alarm;		// stat            ?
    unsigned long
	rtprio,		// stat            real-time priority
	sched,		// stat            scheduling class
	vsize,		// stat            number of pages of virtual memory ...
	rss_rlim,	// stat            resident set size limit?
	flags,		// stat            kernel flags for the process
	min_flt,	// stat            number of minor page faults since process start
	maj_flt,	// stat            number of major page faults since process start
	cmin_flt,	// stat            cumulative min_flt of process and child processes
	cmaj_flt;	// stat            cumulative maj_flt of process and child processes
    char
    	cmd[16];	// stat,status     basename of executable file in call to exec(2)
    int
	pgrp,		// stat            process group id
	session,	// stat            session id
	nlwp,		// stat,status     number of threads, or 0 if no clue
	tty,		// stat            full device number of controlling terminal
        euid, egid,     // stat(),status   effective
	tpgid,		// stat            terminal process group id
	exit_signal,	// stat            might not be SIGCHLD
	processor;      // stat            current (or most recent?) CPU
} proc_stat_t;


// Reads /proc/*/stat files, being careful not to trip over processes with
// names like ":-) 1 2 3 4 5 6".
static void stat2proc(const char* S, proc_stat_t * P) {
    unsigned num;
    char* tmp;

    S = strchr(S, '(') + 1;
    tmp = strrchr(S, ')');
    num = tmp - S;
    if (num >= sizeof P->cmd) 
	num = sizeof P->cmd - 1;
    memcpy(P->cmd, S, num);
    P->cmd[num] = '\0';
    S = tmp + 2;                 // skip ") "

    num = sscanf(S,
		 "%c "
		 "%d %d %d %d %d "
		 "%lu %lu %lu %lu %lu "
		 "%Lu %Lu %Lu %Lu "  /* utime stime cutime cstime */
		 "%ld %ld "
		 "%d "
		 "%ld "
		 "%Lu "  /* start_time */
		 "%lu "
		 "%ld "
		 "%lu %lu %lu %lu %lu %lu "
		 "%*s %*s %*s %*s " /* discard, no RT signals & Linux 2.1 used hex */
		 "%lu %*lu %*lu "
		 "%d %d "
		 "%lu %lu",
		 &P->state,
		 &P->ppid, &P->pgrp, &P->session, &P->tty, &P->tpgid,
		 &P->flags, &P->min_flt, &P->cmin_flt, &P->maj_flt, &P->cmaj_flt,
		 &P->utime, &P->stime, &P->cutime, &P->cstime,
		 &P->priority, &P->nice,
		 &P->nlwp,
		 &P->alarm,
		 &P->start_time,
		 &P->vsize,
		 &P->rss,
		 &P->rss_rlim, &P->start_code, &P->end_code, &P->start_stack, &P->kstk_esp, &P->kstk_eip,
		 /*     P->signal, P->blocked, P->sigignore, P->sigcatch,   */ /* can't use */
		 &P->wchan, /* &P->nswap, &P->cnswap, */  /* nswap and cnswap dead for 2.4.xx and up */
		 /* -- Linux 2.0.35 ends here -- */
		 &P->exit_signal, &P->processor,  /* 2.2.1 ends with "exit_signal" */
		 /* -- Linux 2.2.8 to 2.5.17 end here -- */
		 &P->rtprio, &P->sched  /* both added to 2.5.18 */
			 );

    if(!P->nlwp){
	P->nlwp = 1;
    }
}

static int file2str(const char *directory, const char *what, char *ret, int cap) {
    static char filename[80];
    int fd, num_read;

    sprintf(filename, "%s/%s", directory, what);
    fd = open(filename, O_RDONLY, 0);
    if (fd==-1) 
	return -1;
    num_read = read(fd, ret, cap - 1);
    close(fd);

    if (num_read<=0) 
	return -1;

    ret[num_read] = '\0';
    return num_read;
}

proc_t * get_proc_stats(pid_t pid, proc_t *p) {
    static char path[PATH_MAX], sbuf[1024];
    struct stat statbuf;
    proc_stat_t pf;

    sprintf(path, "/proc/%d", pid);
    if (stat(path, &statbuf)) {
	perror("stat");
	return NULL;
    }

    if (file2str(path, "stat", sbuf, sizeof sbuf) >= 0)
	stat2proc(sbuf, &pf);	/* parse /proc/#/stat */

    p->utime      = pf.utime;
    p->stime      = pf.stime;
    p->start_time = pf.start_time;
    p->rss        = pf.rss*(getpagesize() >> 10);  /* KB */

    return p;
}

