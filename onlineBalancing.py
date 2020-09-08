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

t = 0
time_res = 0.05
t_temp = TIME = 0
steps = 0

(X_l, Y_l, Z_l) = pos_l = pos_l_temp = (0, 0, 0) 
(X_r, Y_r, Z_r) = pos_r = pos_r_temp = (0, 0, 0) 

X_l_traverse, Y_l_traverse, Z_l_traverse = 0, 0, 0
X_r_traverse, Y_r_traverse, Z_r_traverse = 0, 0, 0

foot_left  = (0, 0, 0)
foot_right = (0, 0, 0)

track = []

frxi = foot_left[1]
fryi = 0
frzi = 0

flxi = foot_right[1]
flyi = 0
flzi = 0

th = [0,0,0,0]

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

  return steps, (X_traverse, Y_traverse, Z_traverse)


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


def foot_placement(pos):
  global foot_left,foot_right,pos_l,pos_r,pos_l_temp,pos_r_temp
  if steps%2 == 0:
    pos_r = pos_r_temp 
    foot_left = stepl(pos_l, pos)
    pos_l_temp = pos
  else:
    pos_l = pos_l_temp
    foot_right = stepr(pos_r, pos)
    pos_r_temp = pos
  

def pos():
  global steps#,track
  x = 10*sin(steps+1)
  y = abs(5*sin(steps+1))*1*sign(cos((pi/1)*steps))
  posit = x, y, 2, 0.25
  print(steps)
  #print(posit)
  #track.append( [posit[0], posit[1]] )
  return posit


def step_mark(xl,yl,color):
  global foot_size
  foot_size = 2
  lx = [xl-1.5 ,xl+1.5, xl+1.5, xl-1.5]
  ly = [yl-2 ,yl-2, yl+2, yl+2]
  lz = [0,0,0,0]
  verts = [list(zip(lx,ly,lz))]
  ax.add_collection3d(Poly3DCollection(verts,alpha = 0.3,color = color))


def torqueCalc(theta1,theta2,theta1d,theta2d,theta1dd,theta2dd):
  qdd = mat([ [theta1dd], [theta2dd] ])

  #INERTIAL TERM
  I11 = (m1 + m2)*l_2**2 + m2*(l_2**2) + 2*m2*l_1*l_2*cos(theta2)
  I12 = m2*(l_2**2) + m2*l_1*l_2*cos(theta2) 
  I21 = m2*(l_2**2) + m2*l_1*l_2*cos(theta2) 
  I22 = m2*(l_2**2) 
  I = mat([[ I11, I12 ], [ I21, I22 ]])
  #print I

  #CORIOLIS TERM
  H1 = -m2*l_1*l_2*(2*theta1d*theta2d + theta2d**2)*sin(theta2)
  H2 = -m2*l_2*l_2*theta1d*theta2d*sin(theta2)
  H = mat([ [ H1 ], [ H2 ] ])
  #print H

  #GRAVITY TERM
  G1 = -(m1+m2)*g*l_1*sin(theta1) - m2*g*l_2*sin(theta1+theta2)
  G2 = -m2*g*l_2*sin(theta1 + theta2) 
  G  = mat([ [ G1 ], [ G2 ] ])
  #print G

  #JOINT TORQUES
  Tf = I*qdd + H + G 
  #print Tf
  #print "\nJoint1_torque :",float(Tf[0]),"N/m"
  #print "\nJoint2_torque :",float(Tf[1]),"N/m" 
  return Tf


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


#######################################################################################
lead = 0
#l = 2
while(t < 5):   
    # pos0 = 0, 0, 2, 0.5
    # pos1 = 5*l, 0, 2, 0.5
    # pos2 = 10*l, 0, 2, 0.5
    # pos3 = 15*l, 0, 2, 0.5
    # pos4 = 20*l, 0, 2, 0.5
    # pos5 = 25*l, 0, 2, 0.5
    # pos6 = 30*l, 0, 2, 0.5
    # pos7 = 35*l, 0, 2, 0.5
    # pos8 = 40*l, 0, 2, 0.5
    target = pos()
    foot_placement(target)

    # if steps == 0:
    #   foot_left = stepl(pos0, pos1)
    # if steps == 1:
    #   foot_right = stepr(pos0, pos2)

    # if steps == 2:
    #   foot_left = stepl(pos1, pos3)
    # if steps == 3:
    #   foot_right = stepr(pos2, pos4)

    # if steps == 4:
    #    foot_left = stepl(pos3, pos5)
    # if steps == 5:
    #    foot_right = stepr(pos4, pos6)
    
    # if steps == 6:
    #    foot_left = stepl(pos5, pos7)
    # if steps == 7:
    #    foot_right = stepr(pos6, pos8)
       
    # if steps == 8:
    #    foot_left = stepl(pos7, pos8)

    frx = foot_left[1] 
    fry = foot_left[0]
    frz = foot_left[2]
    
    flx = foot_right[1] 
    fly = foot_right[0]
    flz = foot_right[2]

    #print(flx,fly)

    waist = [0,0,0]

    px =   -2 +flxi - flx +(flx+frx)/2 +waist[0] #+ side_lean*(sin((pi/0.5)*t + pi/2)) -2
    py =   0+fly      -(fly+fry)/2     +waist[1] #-lead*0.2
    pz =   hi-flz                      +waist[2] #- waist_h*abs(cos((pi/0.5)*t))

    #print(pz)
    
    pxr =  2 +frxi - frx +(flx+frx)/2 +waist[0] #+ side_lean*(sin((pi/0.5)*t + pi/2)) +2
    pyr =  0+fry      -(fly+fry)/2    +waist[1] #- lead*0.2
    pzr =  hi-frz                     +waist[2] #- waist_h*abs(cos((pi/0.5)*t))
    print(lead)
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

    ax.plot([ 20] ,[-20 + px] ,[ 0] )
    ax.plot([-20] ,[20 + px] ,[30] )

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
    
    CX,CY,CZ = [(x1+x2)*0.5,(x2+x3)*0.5,(x3+x4)*0.5,(x4+x5)*0.5,(x5+x6+x5r+x6r)*0.25,(x1r+x2r)*0.5,(x2r+x3r)*0.5,(x3r+x4r)*0.5,(x4r+x5r)*0.5], [(y1+y2)*0.5,(y2+y3)*0.5,(y3+y4)*0.5,(y4+y5)*0.5,(y5+y6+y5r+y6r)*0.25,(y1r+y2r)*0.5,(y2r+y3r)*0.5,(y3r+y4r)*0.5,(y4r+y5r)*0.5], [(z1+z2)*0.5,(z2+z3)*0.5,(z3+z4)*0.5,(z4+z5)*0.5,(z5+z6+z5r+z6r)*0.25,(z1r+z2r)*0.5,(z2r+z3r)*0.5,(z3r+z4r)*0.5,(z4r+z5r)*0.5]
    ax.scatter(CX,CY,CZ,c = 'blue',marker = 'o')
    
    CMX,CMY,CMZ = sum(CX)/9,sum(CY)/9,sum(CZ)/9
    ax.plot([CMX],[CMY],[CMZ],c = 'black',marker = 'X',markersize = 5)
    
    zgx = ravel(g*mat(CX)*mat(M).T)/ravel(g*sum(M))
    zgy = ravel(g*mat(CY)*mat(M).T)/ravel(g*sum(M))
    ax.plot([zgx],[zgy],[0],c = 'black',marker = 'X',markersize = 10)

    lead = (target[0]-CMX)

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
   
    step_mark(0,0,'red')
    step_mark(0,0,'green')

    if steps%2 == 0:
      step_mark(pos()[1], pos()[0],'red')
    if steps%2 != 0:
      step_mark(pos()[1], pos()[0],'green')

    bx = [x5, x6, x6r, x5r]
    by = [y5, y6, y6r, y5r]
    bz = [z5, z6, z6r, z5r]
    ver = [list(zip(bx,by,bz))]
    ax.add_collection3d(Poly3DCollection(ver,alpha = 0.3,color = 'blue'))

#===============================================================================
    
    # step_mark(foot_left[1]+5,  foot_left[0],'purple')
    # step_mark(foot_right[1]+5, foot_right[0],'purple')

    # step_mark(foot_left[1],  foot_left[0]+5,'purple')
    # step_mark(foot_right[1], foot_right[0]+5,'purple')

    # step_mark(foot_left[1]-5,  foot_left[0],'purple')
    # step_mark(foot_right[1]-5, foot_right[0],'purple')

    # step_mark(foot_left[1],  foot_left[0]-5,'purple')
    # step_mark(foot_right[1], foot_right[0]-5,'purple')

#===============================================================================
    t+=time_res
    ax.grid()
    plt.pause(0.000000000000000001)
    ax.cla()