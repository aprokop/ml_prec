#!/usr/bin/python
# Plot graph showing residual norms for each level
from numpy import loadtxt
from pylab import *
import matplotlib.pyplot as pyplot

if size(sys.argv) < 2:
    sys.exit('Usage: ./norm.py <trace_file> [output_file]')

norms = loadtxt(sys.argv[1])
n = shape(norms)[0]

step = 0.0
for i in range(1, n):
    if norms[i][0] == 0:
	step = 1./i
	break

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)

colors = { 0:'r', 1:'g', 2:'c', 3:'m', 4:'y', 5:'k', 6:'w' }
# colors = { 0:'1.0', 1:'0.9', 2:'0.8', 3:'0.7', 4:'0.6', 5:'0.5', 6:'0.4' }

ax.plot(arange(0, step*n, step), norms[:,1], '-')
for i in range(0, n):
    ax.plot(i*step, norms[i][1], 'o', color = colors[norms[i][0]])

if 1:
    ax.set_yscale('log')
ax.set_ylim(pow(10,floor(log10(min(norms[:,1])))), pow(10,ceil(log10(max(norms[:,1])))))


if size(sys.argv) < 3:
    show()
else:
    plt.savefig(sys.argv[2])
