#!/usr/bin/python
import numpy
from scipy import interpolate
import sys
from matplotlib.pylab import *

def usage():
    print('image.py <filename> <xy|xz|yz> <coordinate>')

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = linspace(0,1.,N)
    # N+1 indices
    indices = linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

if size(sys.argv) < 4:
    usage()
    sys.exit()

filename = sys.argv[1]
plane = sys.argv[2] # Which 2D plane to output (xy, xz, yz)
coord = sys.argv[3] # Value of the fixed coordinate (z, y, x)

shape = (85,220,60)
# IMPORTANT
# x corresponds to the last dimension
# y corresponds to the second dimension
# z corresponds to the first dimension
# So if we have in C++ x[NODE(i,j,k)] here we have x[k,j,i]

# Set kx and kz
# kx = numpy.fromfile("kx.bin").reshape(shape)
kz = numpy.fromfile("kz.bin").reshape(shape)

# Set x
x = numpy.fromfile(filename)
if False:
    # Clean vector
    # order = -8
    order = -3
    for i in range(len(x)):
	if x[i] < order:
	    x[i] = -100
x = x.reshape(shape)

# color = cm.jet | cm.bone | cm.gray

color = cm.jet
# color = cmap_discretize(cm.jet,7)
if plane == 'xy':
    subplot(121)
    # imshow(kx[coord,:,:], cmap=color)
    imshow(kz[coord,:,:], cmap=color)
elif plane == 'xz':
    subplot(121)
    # imshow(kx[:,coord,:], cmap=color)
    imshow(kz[:,coord,:], cmap=color)
elif plane == 'yz':
    subplot(211)
    # imshow(kx[:,:,coord], cmap=color)
    imshow(kz[:,:,coord], cmap=color)
else:
    print('Unknown value of plane: %s', plane)
    sys.exit(2)
# title('KX: ' + plane + ', fix = ' + coord)
title('KZ: ' + plane + ', fix = ' + coord)
colorbar(ticks=range(10,-10,-1))

color = cmap_discretize(cm.jet,15)
# color = cmap_discretize(cm.gray,15)
if plane == 'xy':
    subplot(122)
    imshow(x[coord,:,:], alpha = 1.0, cmap=color)
elif plane == 'xz':
    subplot(122)
    imshow(x[:,coord,:], alpha = 1.0, cmap=color)
elif plane == 'yz':
    subplot(212)
    imshow(x[:,:,coord], alpha = 1.0, cmap=color)
else:
    print('Unknown value of plane: %s', plane)
    sys.exit(2)

colorbar(ticks=range(0,-100,-1))
title('X: ' + plane + ', fix = ' + coord)

show()
