#!/usr/bin/python
import numpy
import sys
from matplotlib.pylab import *

def usage():
    print('image.py <filename> <xy|xz|yz> <coordinate>')

if size(sys.argv) < 4:
    usage()
    sys.exit()

filename = sys.argv[1]
plane = sys.argv[2] # Which 2D plane to output (xy, xz, yz)
coord = sys.argv[3] # Value of the fixed coordinate (z, y, x)

# Tricky:
# x corresponds to the last dimension
# y corresponds to the second dimension
# z corresponds to the first dimension
# So if we have in C++ x[NODE(i,j,k)] here we have x[k,j,i]
x = numpy.loadtxt(filename)
x = reshape(x, (85,220,60))

# x = reshape(x, (60,220,85), 'f')
# x = swapaxes(x, 0, 2)


if plane == 'xy':
    matshow(x[coord,:,:])
elif plane == 'xz':
    matshow(x[:,coord,:])
elif plane == 'yz':
    matshow(x[:,:,coord])
else:
    print('Unknown value of plane: %s', plane)
    sys.exit(2)

colorbar()
title(plane + ', fix = ' + coord)

show()
