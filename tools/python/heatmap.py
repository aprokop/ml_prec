#!/usr/bin/python

import sys
import numpy
from scipy import interpolate
from matplotlib.pylab import *

from cmap_tools import cmap_discretize

def usage():
    print('Usage: ./heatmap.py <nx> <ny> <filename>')


if size(sys.argv) < 4:
    usage()
    sys.exit()

nx = int(sys.argv[1])
ny = int(sys.argv[2])
filename = sys.argv[3]

shape = (ny, nx)

# Read 1D array
x = numpy.fromfile(filename)
x = x.reshape(shape)

# color = cm.jet | cm.bone | cm.gray
color = cm.gray
# color = cmap_discretize(cm.gray,15)

imshow(x, alpha = 1.0, cmap=color, origin='lower')
# colorbar(ticks=range(-5,0))
colorbar()
title('Heatmap')

show()
