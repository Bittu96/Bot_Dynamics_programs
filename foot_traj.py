import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)

#######################################################################################
t = 0
T = 1

T = 1
H = 5
X = 5
Y = 2

a = (8*H)/(T**3)
b = -(16*H)/(T**2)
c = (8*H)/(T)
d = 0

Z_traverse = a*(t**3) + b*(t**2) + c*t
X_traverse = X*(((t/T) - (1/(2*pi))*cos((2*pi*t)/T - pi/2 )))
Y_traverse = Y*(((t/T) - (1/(2*pi))*cos((2*pi*t)/T - pi/2 )))

speed = 0.4

X_l_traverse, Y_l_traverse, Z_l_traverse = 0, 0, 0
X_r_traverse, Y_r_traverse, Z_r_traverse = 0, 0, 0

track_l = []
track_r = []

def trail(track):
  traillen = 7
  numtrack = len(track)
  if len(track) > traillen:
    numtrack = traillen
  for j in range(0,numtrack):
    ax.scatter(track[len(track)-j-1][0][0],track[len(track)-j-1][1][0],track[len(track)-j-1][2][0],c = 'black',marker = 'o', alpha = 0.5-(j*0.05))

TIME = 0
steps = 0
case = 0

# def forward_step(case):
impact_value = 0.5
T = T/impact_value
Step_length = 2*impact_value

def front(impact):
  steps = 0

while(True):
    steps = int(t/T)
    print("steps :", steps)

    ax.plot([ 30] ,[-30] ,[ -1] )
    ax.plot([-30] ,[ 30] ,[ 20] )

    if steps == 0:
      case = 1
    if steps == 1:
      case = 2
    if steps > 1:
      case = 0

    if case == 0:
      T = 1
      if steps%2 == 0:
        Z_l_traverse = (2*sin((2*pi*t)/T - pi/2 ) + 2)
        Z_r_traverse = 0
      if steps%2 != 0:
        Z_r_traverse = (2*sin((2*pi*t)/T - pi/2 ) + 2)
        Z_l_traverse = 0

    if case == 1:
      X_l_traverse += Step_length*(1+sin((2*pi*t)/T - pi/2 ))*speed*2
      Z_r_traverse = 0
      Z_l_traverse = (2*sin((2*pi*t)/T - pi/2 ) + 2)*impact_value
    
    if case == 2:
      X_r_traverse += Step_length*(1+sin((2*pi*t)/T - pi/2 ))*speed*2
      Z_l_traverse = 0
      Z_r_traverse = (2*sin((2*pi*t)/T - pi/2 ) + 2)*impact_value
    
    if case == -1:
      X_l_traverse -= Step_length*(1+sin((2*pi*t)/T - pi/2 ))
      Z_r_traverse = 0
      Z_l_traverse = (2*sin((2*pi*t)/T - pi/2 ) + 2)*impact_value
    
    if case == -2:
      X_r_traverse -= Step_length*(1+sin((2*pi*t)/T - pi/2 ))
      Z_l_traverse = 0
      Z_r_traverse = (2*sin((2*pi*t)/T - pi/2 ) + 2)*impact_value

    X_l,Y_l,Z_l = [X_l_traverse],[-7.5],[Z_l_traverse]
    X_r,Y_r,Z_r = [X_r_traverse],[ 7.5],[Z_r_traverse]

    track_l.append([X_l,Y_l,Z_l])
    track_r.append([X_r,Y_r,Z_r])

    ax.scatter(X_l,Y_l,Z_l,c = 'red',marker = 'o')
    ax.scatter(X_r,Y_r,Z_r,c = 'green',marker = 'o')
    
    # line, = ax.plot([x6,x6r], [y6,y6r],[z6,z6r], 'black'  , lw=2)
    # line, = ax.plot([x5,x5r], [y5,y5r],[z5,z5r], 'black'  , lw=1)
    #print(track[t-1])
    trail(track_l)
    trail(track_r)

    ax.grid()
    plt.pause(0.000000000000000001)
    ax.cla()

    print(TIME)
    t += speed
    TIME += 0.175