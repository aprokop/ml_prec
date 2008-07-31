#!/bin/bash
if [ $# -lt 1 ]
then
    echo "Usage: ./check_memory.sh <program_name>"
    exit
fi

PROGRAM=check_memory
STATFILE=".stat"
MSGFILE=".msg"

trap "pkill $PROGRAM" SIGINT SIGTERM

rm -f $STATFILE

GSCRIPT=.gscript
echo \
"set terminal postscript size 13,9; 
set output \"stat.ps\";
set grid noxtics ytics x2tics;
set xtics nomirror;
unset x2tics;
set x2tics nomirror rotate;
set x2tics (\"\" 0.001); 
set xlabel \"sec\";
set ylabel \"KB\";" > $GSCRIPT

# execute program; it also may write to $GSCRIPT
`dirname $0`/check_memory/check_memory $*

echo \
"set style line 1 lt rgb \"red\" lw 1;
set style line 2 lt rgb \"black\" lw 2;
plot \"$STATFILE\" with filledcurves below ls 1 title \"memory\", \"$STATFILE\" with lines ls 2" >> $GSCRIPT

gnuplot $GSCRIPT &>/dev/null

rm -f $GSCRIPT
rm -f $STATFILE
rm -f $MSGFILE
