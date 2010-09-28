#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

N = 13
# matrix-3D-869K-cfl1-alpha1.0
# values = ( 0, 4872, 56778, 338711, 403540, 63353, 600, 255, 229, 322, 564, 42, 43 )
# matrix-3D-869K-cfl1-alpha0.1
# values = ( 0, 0, 4872, 56889, 339670, 408888, 51638, 3316, 885, 938, 1477, 680, 56 )
# matrix-3D-869K-cfl1-alpha0.01
# values = ( 0, 0, 0, 4875, 58088, 350799, 161344, 115838, 74297, 37966, 63684, 2362, 56 )
# matrix-3D-869K-cfl1-alpha0.001
# values = ( 0, 0, 0, 0, 4905, 70432, 97657, 100262, 75372, 51608, 395231, 73786, 56 )
# spe, unsymmetric, shift=0.1
values = ( 0, 0, 204, 44537, 266367, 120032, 58991, 45252, 40568, 40613, 281932, 223504, 0 )


ticks  = ('1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '0.2', '0.3', '0.4', '0.5', '0.9', '1', 1.0000001)

ind = np.arange(N)  # the x locations for the groups
width = 1         # the width of the bars

rects = plt.bar(ind, values, width, color='r')

# add some
plt.ylabel('')
plt.title('')
plt.xticks(ind+width, ticks )

plt.title('matrix-3D-869K')

def autolabel(rects):
    # attach some text labels
    for rect in rects:
	height = rect.get_height()
	plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
		 ha='center', va='bottom')

autolabel(rects)

plt.show()
