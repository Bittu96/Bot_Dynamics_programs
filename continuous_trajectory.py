import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import time

fig = plt.figure()
ax = Axes3D(fig)

#######################################################################################

t = 0
time_res = 0.05
t_temp = TIME = 0

steps = 0

(X_l, Y_l, Z_l) = (0, 0, 0)
foot_left_i  = [ 3, 0, 0]


def step_mark(xl,yl,color):
  global foot_size
  foot_size = 2
  lx = [xl-1.5 ,xl+1.5, xl+1.5, xl-1.5]
  ly = [yl-2 ,yl-2, yl+2, yl+2]
  lz = [0,0,0,0]

  verts = [list(zip(lx,ly,lz))]
  ax.add_collection3d(Poly3DCollection(verts,alpha = 0.3,color = color))


def Trajectory(X,Y,H,T):
  global t_temp,track_l,trail,steps,X_traverse,Y_traverse,Z_traverse
  a = (8*H)/(T**3)
  b = -(16*H)/(T**2)
  c = (8*H)/(T)
  d = 0
  if t_temp <= T:
    Z_traverse = a*(t_temp**3) + b*(t_temp**2) + c*t_temp
    X_traverse = X*(((t_temp/T) - (1/(2*pi))*cos((2*pi*t_temp)/T - pi/2 )))
    Y_traverse = Y*(((t_temp/T) - (1/(2*pi))*cos((2*pi*t_temp)/T - pi/2 )))
    # ax.scatter(X_traverse,Y_traverse,Z_traverse,c = 'red',marker = 'o')
    t_temp+=time_res
    # track_l.append([[X_traverse],[Y_traverse],[Z_traverse]])
    # print(track_l)
    #trail(track_l,'red')
  else:
    t_temp = 0
    steps += 1
    
  return steps, [X_traverse, Y_traverse, Z_traverse]


def step(pos0,pos,leg):
  global X_l, Y_l, Z_l, X_r, Y_r, Z_r
  if steps >= 0:
    for i in range(0,1):
      #Trajectory(5,5,2,1)
      traject = Trajectory(pos[0]-pos0[0] ,pos[1]-pos0[1], pos[2], pos[3])[1]
      #print(traject)
  
      if leg == 'L':
        X_l = traject[0]
        Y_l = traject[1]
        Z_l = traject[2]

  return [X_l + pos0[0], Y_l + pos0[1], Z_l]

while(t < 2.3):
    ax.plot([ 30] ,[-30] ,[ -1] )
    ax.plot([-30] ,[ 30] ,[ 20] )

    pos0 = 0 ,0  ,2 ,0.5
    pos1 = 0 ,10 ,2 ,0.5
    pos2 = 10,10 ,2 ,0.5
    pos3 = 10,0  ,2 ,0.5

    print(sign(0))
    if steps == 0:
      ff = step(pos0, pos1, 'L')
    if steps == 1:
      ff = step(pos1, pos2, 'L')
    if steps == 2:
      ff = step(pos2, pos3, 'L')
    if steps == 3:
      ff = step(pos3, pos0, 'L')
    print(ff)
    foot_left  = ff

    print("steps : ", steps)
    ax.scatter(foot_left[0] , foot_left[1], foot_left[2], c = 'red',marker = 'o')
    step_mark(0,0,'green')
    step_mark(0,10,'green')
    step_mark(10,0,'green')
    step_mark(10,10,'green')
    ax.grid()
    plt.pause(0.000000000000000001)
    ax.cla()

    t += time_res
    TIME += 0.175