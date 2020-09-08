import time
import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)

def T(p,q,r,x,y,z):
    Rz = mat([ [cos(r), -sin(r),       0], [sin(r), cos(r),      0], [0      , 0     ,      1]])
    Ry = mat([ [cos(q), 0      ,  sin(q)], [0     ,      1,      0], [-sin(q), 0     , cos(q)]])
    Rx = mat([ [1     , 0      ,       0], [0     , cos(p),-sin(p)], [0      ,sin(p) , cos(p)]])
    #Rxyz = Rx*Ry*Rz
    Rxyz = Rz*Ry*Rx
    temp = ravel(Rxyz).T
    Ti = mat([ 
    [temp[0],temp[1],temp[2],x],
    [temp[3],temp[4],temp[5],y],
    [temp[6],temp[7],temp[8],z],
    [0      ,0      ,0      ,1] 
    ])
    return Ti

def iksol(X,Y,Z):
  theta1 = atan2(X,Z)

  _T2 = T(0,0,0,0,0,0) * T(0,theta1-pi/2,0,0,0,0) * T(0,0,0,a1,0,0)
  x2 = _T2[0,3]
  y2 = _T2[1,3]
  z2 = _T2[2,3]

  l = sqrt( (x2-X)**2 + (y2-Y)**2 + (z2-Z)**2 ) 
  lz= sqrt( (x2-X)**2 + (z2-Z)**2 ) 

  c21 = atan2(lz,l) 
  s21 = sqrt(1 - c21**2)
  theta21 = atan2(Y,lz)
  c22 = (a2**2 + l**2 - a3**2)/(2*a2*l) 
  s22 = sqrt(1 - c22**2)
  theta22 = atan2(s22,c22)
  theta2 = theta21 - theta22

  s23 = (a2*s22)/a3
  c23 = sqrt(1 - s23**2)
  theta3 = theta22 + atan2(s23,c23)
  return theta1,theta2,theta3

a1 = 2
a2 = 5
a3 = 5

X = 11
Y = 3
Z = 5

SimulationRunning = True

theta1i = 0#pi/2 + pi/3 #radians
theta2i = 0#0 + pi/3

theta1di = 0 #radians/sec
theta2di = 0

theta1ddi = 0 #radians/sec2
theta2ddi = 0

#\\\\\\\\\\\\\  PARAMETERS \\\\\\\\\\\\\\\\
m1 = 3      #kg
m2 = 3      #kg
l1 = 0.3    #m
l2 = 0.3    #m
g  = 9.81   #m/s2
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

def torqueCalc(theta1,theta2,theta1d,theta2d,theta1dd,theta2dd):
  qdd = mat([ [theta1dd], [theta2dd] ])

  #INERTIAL TERM
  I11 = (m1 + m2)*l2**2 + m2*(l2**2) + 2*m2*l1*l2*cos(theta2)
  I12 = m2*(l2**2) + m2*l1*l2*cos(theta2) 
  I21 = m2*(l2**2) + m2*l1*l2*cos(theta2) 
  I22 = m2*(l2**2) 
  I = mat([[ I11, I12 ], [ I21, I22 ]])
  #print I

  #CORIOLIS TERM
  H1 = -m2*l1*l2*(2*theta1d*theta2d + theta2d**2)*sin(theta2)
  H2 = -m2*l2*l2*theta1d*theta2d*sin(theta2)
  H = mat([ [ H1 ], [ H2 ] ])
  #print H

  #GRAVITY TERM
  G1 = -(m1+m2)*g*l1*sin(theta1) - m2*g*l2*sin(theta1+theta2)
  G2 = -m2*g*l2*sin(theta1 + theta2) 
  G =mat([ [ G1 ], [ G2 ] ])
  #print G

  #JOINT TORQUES
  T = I*qdd + H + G 
  print "\nJoint1_torque :",int(T[0]),"N/m"
  print "\nJoint2_torque :",int(T[1]),"N/m" 

    


####################################### VISUAL #############################################################################
####################################### VISUAL #############################################################################
####################################### VISUAL #############################################################################

time = 0
while SimulationRunning:
  X = 0#- 5*sin(time)
  Y = 0#+ 5*sin(time)
  Z = -10#+ 5*sin(time)

  theta1 = iksol(X,Y,Z)[0]
  theta2 = iksol(X,Y,Z)[1]
  theta3 = iksol(X,Y,Z)[2]

  theta4 = 0

  theta1d = -theta1i + theta1
  theta2d = -theta2i + theta2

  theta1i = theta1
  theta2i = theta2

  theta1dd = -theta1di + theta1d
  theta2dd = -theta2di + theta2d

  theta1di = theta1d
  theta2di = theta2d

  torqueCalc(theta1,theta2,theta1d,theta2d,theta1dd,theta2dd)

  _0T1 = T(0,0,0,0,0,0) * T(0,theta1-pi/2,0,0,0,0)
  _1T2 = T(0,0,theta2,a1,0,0)
  _2T3 = T(0,0,theta3,a2,0,0)
  _3T4 = T(0,0,theta4,a3,0,0)

  _0T1 = _0T1
  _0T2 = _0T1 * _1T2
  _0T3 = _0T2 * _2T3
  _0T4 = _0T3 * _3T4

  x1 = _0T1[0,3]
  y1 = _0T1[1,3]
  z1 = _0T1[2,3]

  x2 = _0T2[0,3]
  y2 = _0T2[1,3]
  z2 = _0T2[2,3]

  x3 = _0T3[0,3]
  y3 = _0T3[1,3]
  z3 = _0T3[2,3]

  x4 = _0T4[0,3]
  y4 = _0T4[1,3]
  z4 = _0T4[2,3]

  ax.plot([X] ,[Y] ,[Z] ,'gs')
  ax.plot([15] ,[-10] ,[0] )
  ax.plot([-10] ,[15] ,[-10] )

  X,Y,Z = [x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4]
  ax.scatter(X,Y,Z,c = 'red',marker = 'o')

  line, = ax.plot([x1,x2], [y1,y2],[z1,z2], 'blue', lw=1)
  line, = ax.plot([x2,x3], [y2,y3],[z2,z3], 'green', lw=1)
  line, = ax.plot([x3,x4], [y3,y4],[z3,z4], 'red', lw=1)

  ax.grid()
  plt.pause(0.000000001)
  ax.cla()
  time += 0.1 