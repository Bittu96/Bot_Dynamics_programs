import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)

#\\\\\\\\\\\\\  PARAMETERS \\\\\\\\\\\\\\\\

a1 = 2
a2 = 9
a3 = 9
a4 = 2
a5 = 10

hi = 18

waist_h = 5
side_lean = 1

m1,m2,m3,m4,m5 = 2,9,9,2,10
M = [m1,m2,m3,m4,m5,m1,m2,m3,m4]
g  = 9.81   #m/s2

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

t= 0
jt= jz =0
time_res = 0.05

foot_left  = (0, 1, 0)
foot_right = (0, -1, 0)

frxi = foot_left[1]
fryi = 0
frzi = 0

flxi = foot_right[1]
flyi = 0
flzi = 0

th = [0,0,0,0]

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
  global T
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

def jump(thrust):
  global jt,jz; c1= c2= 5; lamda1= -3; lamda2= -3
  if jz >= 0:
    #jz= thrust*jt - g*jt**2 
    jz = 5*sin(t)# *(c1*(e**(lamda1*jt)) + c2*(e**(lamda2*jt)) ) 

    jt+= time_res
  else:
    jz = 0
    jt = 0
   
  return jz


#######################################################################################
while(t < 5):   
    frx = foot_left[1] 
    fry = foot_left[0]
    frz = foot_left[2] + jump(10)
    
    flx = foot_right[1] 
    fly = foot_right[0]
    flz = frz#foot_right[2]

    #print(flx,fly)

    waist = [0,0,0]

    px =   flxi - flx +(flx+frx)/2 +waist[0] #+ side_lean*(sin((pi/0.5)*t + pi/2)) -2
    py =   0+fly      -(fly+fry)/2     +waist[1] #-lead*0.2
    pz =   hi -flz                     +waist[2] #- waist_h*abs(cos((pi/0.5)*t))

    #print(pz)
    
    pxr =  frxi - frx +(flx+frx)/2 +waist[0] #+ side_lean*(sin((pi/0.5)*t + pi/2)) +2
    pyr =  0+fry      -(fly+fry)/2    +waist[1] #- lead*0.2
    pzr =  hi -frz                    +waist[2] #- waist_h*abs(cos((pi/0.5)*t))
    
    legl = iksol(px, py, pz)
    legr = iksol(pxr, pyr, pzr)

    theta1 = legl[0]
    theta2 = -legl[1]
    theta3 = -legl[2]
    theta4 = legl[1]+legl[2]
    theta5 = -legl[0]

    theta1r = legr[0]
    theta2r = -legr[1]
    theta3r = -legr[2]
    theta4r = legr[1]+legr[2]
    theta5r = -legr[0]
  
    av1 = theta1 - th[0]
    aa1 = th[0]  - th[1]    
    
    av2 = theta2 - th[2]
    aa2 = th[2]  - th[3]
    #print(theta3*180/pi, theta4*180/pi, theta5*180/pi)
  
    #---------------------------------------Kinematics q-------------------------------------------------------------------
    
    _0T0 = T(0,0,0     ,flx,fly,flz)
    _0T1 = T(0,theta1-pi/2,0,0,0,0)
    _1T2 = T(0,0,theta2,a1,0,0)
    _2T3 = T(0,0,theta3,a2,0,0)
    _3T4 = T(0,0,theta4,a3,0,0)
    _4T5 = T(0,theta5,0,a4,0,0)
    _5T6 = T(0,0,1     ,a5,0,0)

    _0T0_ = T(0,0,0     ,frx,fry,frz)
    _0T1_ = T(0,theta1r-pi/2,0,0,0,0)
    _1T2_ = T(0,0,theta2r,a1,0,0)
    _2T3_ = T(0,0,theta3r,a2,0,0)
    _3T4_ = T(0,0,theta4r,a3,0,0)
    _4T5_ = T(0,theta5r,0,a4,0,0)
    _5T6_ = T(0,0,0     ,a5,0,0)
    
    #//////////////////////////////////////////// D-H Parameters -_---------------------------------------------------

    _0T1 = _0T0 * _0T1
    _0T2 = _0T1 * _1T2
    _0T3 = _0T2 * _2T3
    _0T4 = _0T3 * _3T4
    _0T5 = _0T4 * _4T5
    _0T6 = _0T5 * _5T6

    _0T1_ = _0T0_ * _0T1_
    _0T2_ = _0T1_ * _1T2_
    _0T3_ = _0T2_ * _2T3_
    _0T4_ = _0T3_ * _3T4_
    _0T5_ = _0T4_ * _4T5_
    _0T6_ = _0T5_ * _5T6_

    ax.plot([ 20] ,[-20 + px] ,[ -5] )
    ax.plot([-20] ,[20 + px] ,[50] )

    x1 = _0T1[0,3]
    y1 = _0T1[1,3]
    z1 = _0T1[2,3]
    #axesi.append([x1,y1,z1])
    x2 = _0T2[0,3]
    y2 = _0T2[1,3]
    z2 = _0T2[2,3]
    #axesi.append([x2,y2,z2])
    x3 = _0T3[0,3]
    y3 = _0T3[1,3]
    z3 = _0T3[2,3]
    #axesi.append([x3,y3,z3])
    x4 = _0T4[0,3]
    y4 = _0T4[1,3]
    z4 = _0T4[2,3]
    #axesi.append([x4,y4,z4])
    x5 = _0T5[0,3]
    y5 = _0T5[1,3]
    z5 = _0T5[2,3]

    x6 = _0T6[0,3]
    y6 = _0T6[1,3]
    z6 = _0T6[2,3]

    x1r = _0T1_[0,3]
    y1r = _0T1_[1,3]
    z1r = _0T1_[2,3]
    #axesi.append([x1,y1,z1])
    x2r = _0T2_[0,3]
    y2r = _0T2_[1,3]
    z2r = _0T2_[2,3]
    #axesi.append([x2,y2,z2])
    x3r = _0T3_[0,3]
    y3r = _0T3_[1,3]
    z3r = _0T3_[2,3]
    #axesi.append([x3,y3,z3])
    x4r = _0T4_[0,3]
    y4r = _0T4_[1,3]
    z4r = _0T4_[2,3]
    #axesi.append([x4,y4,z4])
    x5r = _0T5_[0,3]
    y5r = _0T5_[1,3]
    z5r = _0T5_[2,3]

    x6r = _0T6_[0,3]
    y6r = _0T6_[1,3]
    z6r = _0T6_[2,3]


    X,Y,Z = [x1,x2,x3,x4,x5,x6,x1r,x2r,x3r,x4r,x5r,x6r], [y1,y2,y3,y4,y5,y6,y1r,y2r,y3r,y4r,y5r,y6r], [z1,z2,z3,z4,z5,z6,z1r,z2r,z3r,z4r,z5r,z6r]
    ax.scatter(X,Y,Z,c = 'black',marker = 'o')
    
    line, = ax.plot([x1,x2], [y1,y2],[z1,z2], 'black', lw=1)
    line, = ax.plot([x2,x3], [y2,y3],[z2,z3], 'black', lw=1)
    line, = ax.plot([x3,x4], [y3,y4],[z3,z4], 'black', lw=2)
    line, = ax.plot([x4,x5], [y4,y5],[z4,z5], 'black', lw=1)
    line, = ax.plot([x5,x6], [y5,y6],[z5,z6], 'black', lw=1)

    line, = ax.plot([x1r,x2r], [y1r,y2r],[z1r,z2r], 'black', lw=1)
    line, = ax.plot([x2r,x3r], [y2r,y3r],[z2r,z3r], 'black', lw=1)
    line, = ax.plot([x3r,x4r], [y3r,y4r],[z3r,z4r], 'black', lw=2)
    line, = ax.plot([x4r,x5r], [y4r,y5r],[z4r,z5r], 'black', lw=1)
    line, = ax.plot([x5r,x6r], [y5r,y6r],[z5r,z6r], 'black', lw=1)

    line, = ax.plot([x6,x6r], [y6,y6r],[z6,z6r], 'black', lw=2)
    line, = ax.plot([x5,x5r], [y5,y5r],[z5,z5r], 'black', lw=1)

    t+=time_res
    ax.grid()
    plt.pause(0.000000000000000001)
    ax.cla()