#!/usr/bin/python
# Plot graph showing residual norms for each level
from numpy import loadtxt
from pylab import *
import matplotlib.pyplot as pyplot

if size(sys.argv) < 1:
    sys.exit('Usage: ./2level [output_file]')

real = loadtxt("2level_real.dat")
form = loadtxt("2level_form.dat")

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)

ax.plot(real[:,0], real[:,1])
ax.plot(form[:,0], form[:,1])

if size(sys.argv) < 2:
    show()
else:
    plt.savefig(sys.argv[2])
