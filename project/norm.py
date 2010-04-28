#!/usr/bin/python
from numpy import loadtxt
from pylab import *
import matplotlib.pyplot as pyplot

norms = loadtxt(sys.argv[1])

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)

colors = { 0:'r', 1:'g', 2:'c', 3:'m', 4:'y', 5:'k', 6:'w' }
# colors = { 0:'1.0', 1:'0.9', 2:'0.8', 3:'0.7', 4:'0.6', 5:'0.5', 6:'0.4' }

ax.plot(norms[:,1], '-')
for i in range(0, shape(norms)[0]):
    ax.plot(i, norms[i][1], 'o', color = colors[norms[i][0]])

ax.set_yscale('log')

show()
