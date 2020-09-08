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
(X_r, Y_r, Z_r) = (0, 0, 0)

foot_left_i  = [ 3, 0, 0]
foot_right_i = [-3, 0, 0]

track_l = []
track_r = []

def trail(track,color):
  traillen = 3
  numtrack = len(track)
  if len(track) > traillen:
    numtrack = traillen
  for j in range(0,numtrack):
    ax.scatter(track[len(track)-j-1][0][0], track[len(track)-j-1][1][0], track[len(track)-j-1][2][0], c = color, marker = 'o', alpha = 0.5-(j*0.05))


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
    t_temp+=time_res
   
  else:
    t_temp = 0
    steps += 1
    
  return steps, [X_traverse, Y_traverse, Z_traverse]


def stepl(pos0,pos):
  global X_l, Y_l, Z_l, X_r, Y_r, Z_r
  traject = Trajectory(pos[0]-pos0[0] ,pos[1]-pos0[1], pos[2], pos[3])[1]
  X_l = traject[0]
  Y_l = traject[1]
  Z_l = traject[2]
  return [X_l+ pos0[0], Y_l+ pos0[1], Z_l]

def stepr(pos0,pos):
  global X_l, Y_l, Z_l, X_r, Y_r, Z_r
  traject = Trajectory(pos[0]-pos0[0] ,pos[1]-pos0[1], pos[2], pos[3])[1]
  X_r = traject[0]
  Y_r = traject[1]
  Z_r = traject[2]
  return [X_r+ pos0[0], Y_r+ pos0[1], Z_r]

  
def step_mark(xl,yl,color):
  global foot_size
  foot_size = 2
  lx = [xl-1.5 ,xl+1.5, xl+1.5, xl-1.5]
  ly = [yl-2 ,yl-2, yl+2, yl+2]
  lz = [0,0,0,0]

  verts = [list(zip(lx,ly,lz))]
  ax.add_collection3d(Poly3DCollection(verts,alpha = 0.3,color = color))


while(t < 4.5):
    ax.plot([ 30] ,[-30] ,[ -1] )
    ax.plot([-30] ,[ 30] ,[ 20] )

    pos0 = 0 ,0  ,2 ,0.5
    pos1 = 0 ,10 ,2 ,0.5
    pos2 = 10,10 ,2 ,0.5
    pos3 = 10,0  ,2 ,0.5

    if steps == 0:
      ff = stepl(pos0, pos1)
      foot_left  = ff
      foot_right = [X_r, Y_r, Z_r]
    if steps == 1:
      ff = stepr(pos0, pos1)
      foot_right = ff

    if steps == 2:
      ff = stepl(pos1, pos2)
      foot_left  = ff
    if steps == 3:
      ff = stepr(pos1, pos2)
      foot_right = ff

    if steps == 4:
       ff = stepl(pos2, pos3)
       foot_left  = ff
    if steps == 5:
       ff = stepr(pos2, pos3)
       foot_right = ff
    
    if steps == 6:
       ff = stepl(pos3, pos0)
       foot_left  = ff
    if steps == 7:
       ff = stepr(pos3, pos0)
       foot_right = ff
       
    step_mark(0,0,'green')
    step_mark(0,10,'green')
    step_mark(10,0,'green')
    step_mark(10,10,'green')
 
    print("steps : ", steps)
    ax.scatter(foot_left[0] , foot_left[1], foot_left[2], c = 'red',marker = 'o')
    ax.scatter(foot_right[0], foot_right[1], foot_right[2], c = 'green',marker = 'o')
    
    #trail(track_l, 'red')
    #trail(track_r)

    ax.grid()
    plt.pause(0.000000000000000001)
    ax.cla()

    #print(t, TIME)
    t += time_res
    TIME += 0.175