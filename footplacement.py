import time
import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)

SimulationRunning = True
res = 1

ix = 5
iy = 5
iz = 0

sx = -5
sy = -5
sz = 0
left = True

while SimulationRunning:
  x = 0
  y = 0
  if left == True:
    ts = 0
    H = 3
    T = 2
    a1 =  4*H / T
    a2 = -4*H / T**2
    a3 = (sx-ix)/T
    a4 = (sy-iy)/T
    while ts < T:
      x = ix + a3*ts
      y = iy + a4*ts
      z = iz + a1*ts + a2*(ts**2)
      #print z
      ax.plot([x],[y],[z],'ro',markersize = 2)

      ax.plot([ix],[iy],[iz],'bo',markersize = 2)
      ax.plot([sx],[sy],[sz],'go',markersize = 2)
      ax.plot([10] ,[-10] ,[10] )
      ax.plot([-10] ,[10] ,[0] )
      ax.grid()
      plt.pause(0.00000000000001)
      ts += 0.05
  ax.cla()
