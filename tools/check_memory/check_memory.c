#include "proc.h"
#include <assert.h>
#include <signal.h>
#include <string.h>

float time     = 0;
long rss       = 0;
FILE * fscript = NULL;
FILE * fdstat  = NULL;

useconds_t usec = 5000;

static void signal_handler(int signo){
    if (signo == SIGUSR1) {
	FILE * fmsg = fopen(".msg", "r");
	char * msg = NULL;
	size_t len = 0;

	if (fmsg != NULL) {
	    if (getline(&msg, &len, fmsg) != -1) {
		msg[strlen(msg)-1] = 0;
		fprintf(fscript, "set x2tics add (\"%s (%.2f, %ld)\" %.3f);\n", msg, time, rss, time);
		fflush(fscript);
		if (msg)
		    free(msg);
	    } else {
		perror("getline");
	    }
	    fclose(fmsg);

	} else {
	    perror("fopen");
	    fprintf(stderr, "Caught SIGUSR1 but couldn't open \".msg\"\n");
	}
    }
}

int main(int argc, char *argv[]) {
    char **newargv = argv+1;
    pid_t cpid, w;
    int status;
    proc_t pstat;
    struct sigaction sa;

    if (argc < 3) {
	printf("Usage: %s <max_memory_size_in_MB> <program_name> [<program_params>]\n", argv[0]);
	exit(EXIT_SUCCESS);
    }
    uint max_mem_mb = atoi(argv[1]);

    /* set signal handlers */
    do {
	struct sigaction sa;
	int i = 32;
	memset(&sa, 0, sizeof(sa));
	sa.sa_handler = signal_handler;
	sigfillset(&sa.sa_mask);
	while(i--) switch(i){
	    default:
		sigaction(i,&sa,NULL);
	    case 0:
	    case SIGINT:   /* ^C */
	    case SIGTSTP:  /* ^Z */
	    case SIGTTOU:  /* see stty(1) man page */
	    case SIGQUIT:  /* ^\ */
	    case SIGPROF:  /* profiling */
	    case SIGKILL:  /* can not catch */
	    case SIGSTOP:  /* can not catch */
	    case SIGWINCH: /* don't care if window size changes */
		;
	}
    } while (0);

    fdstat = fopen(".stat", "w");
    fscript = fopen(".gscript", "a");
    if (fdstat == NULL || fscript == NULL) {
	perror("fopen");
	exit(EXIT_FAILURE);
    }
    fprintf(fdstat, "0 0\n");

    /* fork process */
    cpid = fork();
    if (cpid == -1) {
	perror("fork");
	exit(EXIT_FAILURE);
    }

    if (cpid == 0) {            /* Code executed by child */
	char *newenviron[] = { NULL };

	// printf("Child PID is %ld\n", (long) getpid());
	execve(argv[2], newargv, newenviron);

	perror("execve");   /* execve() only returns on error */
	exit(EXIT_FAILURE);

    } else {                    /* Code executed by parent */
	while (1) {
	    /* sleep */
	    usleep(usec);

	    w = waitpid(cpid, &status, WNOHANG);
	    if (w == -1) {
		perror("waitpid");
		exit(EXIT_FAILURE);
	    }

	    time += usec*1e-6;

	    if (w > 0) {
		fprintf(fdstat, "%.3f 0\n", time);
		break;
	    }

	    if (get_proc_stats(cpid, &pstat) == NULL)
		assert(0);

	    rss = pstat.rss;
	    fprintf(fdstat, "%.3f %ld\n", time, rss);

	    if (rss > max_mem_mb*1024) {
		kill(cpid, SIGKILL);
		fprintf(stderr, "Maximum memory size (%dMB) exceeded, terminating...\n", max_mem_mb);
	    }
	};
    }
    fclose(fdstat);

    fprintf(fscript, "set xrange  [0:%.3f];\n", time);
    fprintf(fscript, "set x2range [0:%.3f];\n", time);
    fclose(fscript);

    exit(EXIT_SUCCESS);
}

