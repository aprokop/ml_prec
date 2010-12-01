#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

N = 13
# values = (...)
values = ( 0, 0, 204, 44537, 266367, 120032, 58991, 45252, 40568, 40613, 281932, 223504, 0 )

ticks  = ('1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '0.2', '0.3', '0.4', '0.5', '0.9', '1', 1.0000001)

ind = np.arange(N)  # the x locations for the groups
width = 1         # the width of the bars

rects = plt.bar(ind, values, width, color='r')

# add some
plt.ylabel('')
plt.title('')
plt.xticks(ind+width, ticks )

plt.title('spe')

def autolabel(rects):
    # attach some text labels
    for rect in rects:
	height = rect.get_height()
	plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
		 ha='center', va='bottom')

autolabel(rects)

plt.show()
