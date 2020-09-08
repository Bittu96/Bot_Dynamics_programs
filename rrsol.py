import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)

H = 8
L = 3
d1 = d2 = 10

A = ((H**2 + L**2 + d1**2 - d2**2)/(2*d1*sqrt(H**2 + L**2)))
B = sqrt(1 - A**2)

t1 = atan2(A,B) - atan2(H,L)
t2 = atan2(L - d1*sin(t1), H - d2*cos(t1))
t3 = 0

X = d1*cos(t1) + d2*cos(t2)
Y = d1*sin(t1) + d2*sin(t2)
Z = 0

x = d1*cos(t1)
y = d1*sin(t1)
z = 0#d1*cos(t3)

#---------------------------------------plotting section-------------------------------------------------------------------
ax.plot([10] ,[10] ,[5] )
ax.plot([10],[10],[-20])
ax.plot([x],[y],[z],'ro')
ax.plot([X],[Y],[Z],'ro')

line, = ax.plot([0,x], [0,y],[0,z], 'black', lw=2)
line, = ax.plot([x,X], [y,Y],[z,Z], 'black', lw=2)

plt.grid()
plt.show()
