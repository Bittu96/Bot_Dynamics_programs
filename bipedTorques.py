import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)
#ax  = fig.add_subplot(2, 1, 1, projection='3d')
#ax1 = fig.add_subplot(2, 1, 2)
#ax2 = fig.add_subplot(2, 2, 1)
#ax3 = fig.add_subplot(2, 2, 2)

fig, axes = plt.subplots(nrows=2, ncols=2)
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

firstStep = True
secondStepLoop = True

a1 = 1
a2 = 9
a3 = 8
a4 = 1
a5 = 5

#\\\\\\\\\\\\\  PARAMETERS \\\\\\\\\\\\\\\\
m1 = 0.9      #kg
m2 = 0.8      #kg
l_1 = 0.9    #m
l_2 = 0.8    #m
g  = 9.81   #m/s2
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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
  print Tf
  print "\nJoint1_torque :",float(Tf[0]),"N/m"
  print "\nJoint2_torque :",float(Tf[1]),"N/m" 
  return Tf

def iksol(X,Y,Z):
  theta1 = atan2(X,Z)
  global T
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

iksol(0,0,0)
#######################################################################################
t = 0
y = 0
x = 0

T = 0.5
H = 0.7

a = -(8*H)/(T**3)
b =  (8*H)/(T**2)

c = 0
d = 0

i=0

l1 = 0
l2 = 0
frx = 1
fry = 0
frz = 0

flx = -1
fly = 0
flz = 0

hi = 15
yi = 2
speed = 1
th = [0,0,0,0]

while(firstStep):		 
	i += 0.04*speed
	t += 0.04*speed
	
	if i<1:
		x = 5*i
		l = 1
	else:
		x = 5*i
		l = 2
		
	#print i
	y1 = (sin(i*2*pi) ) * ( sign( sin(i*2*pi*0.5) ) + sign( sin(i*2*pi*0.5 +pi/2) ) + sign( sin(i*2*pi) ) + 1)*0.5
	y2 = (sin(i*2*pi) ) * ((1 + sign(sin(i*2*pi)) + sign( sin(i*2*pi*0.5 +pi) ) + sign( sin(i*2*pi*0.5 -pi/2) ) ))*0.5

	l1 += l*0.1*speed*( sign( sin(i*pi) )      + sign( sin(i*pi +   pi/2) ) + sign( sin(i*2*pi) ) + 1 )
	l2 += l*0.1*speed*( sign( sin(i*pi + pi) ) + sign( sin(i*pi + 1.5*pi) ) + sign( sin(i*2*pi + 2*pi) ) + 1 )
    
	frx = 1
	fry = 0+l2
	frz = 0+y2

	flx = -1
	fly =  0+l1
	flz =  0+y1

	px =    0-flx-1 + 0.5*cos(i*pi - pi/4)
	py =   yi+fly -x#- 5*i + 5*(1/pi)*sin(i*pi)
	pz =   hi-flz
	
	pxr =   0-frx+1 + 0.5*cos(i*pi - pi/4)
	pyr =  yi+fry -x#- 5*i + 5*(1/pi)*sin(i*pi)
	pzr =  hi-frz
	
	legl = iksol(px, py, pz)
	legr = iksol(pxr, pyr, pzr)

	theta1 = legl[0]
	theta2 = -legl[1]-speed*0.01
	theta3 = -legl[2]
	theta4 = legl[1]+legl[2] + pi/6*(1+sin(pi*i*2))*0.25/speed 
	theta5 = -legl[0]

	theta1r = legr[0]
	theta2r = -legr[1]-speed*0.01
	theta3r = -legr[2]
	theta4r = legr[1]+legr[2] + pi/6*(1+sin(pi*i*2))*0.25/speed
	theta5r = -legr[0]
	
	print theta1,theta2,theta3,theta4,theta5
	av1 = theta1 - th[0]
	aa1 = th[0]  - th[1]
	
	av2 = theta2 - th[2]
	aa2 = th[2]  - th[3]
	#print theta1,theta2, th[0],th[2], th[1],th[3]
	Tq = torqueCalc(theta1,theta2, av1, av2, aa1, aa2)
	
	#---------------------------------------Kinematics section-------------------------------------------------------------------

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
	
	_0T0 = T(0,0,0     ,flx,fly,flz)
	_0T1 = T(0,theta1-pi/2,0,0,0,0)
	_1T2 = T(0,0,theta2,a1,0,0)
	_2T3 = T(0,0,theta3,a2,0,0)
	_3T4 = T(0,0,theta4,a3,0,0)
	_4T5 = T(0,theta5,0,a4,0,0)
	_5T6 = T(0,0,0     ,a5,0,0)

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



	ax.plot([10] ,[-5] ,[0] )
	ax.plot([-10] ,[80] ,[20] )

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

	X,Y,Z = [0,x1,x2,x3,x4,x5,x6,x1r,x2r,x3r,x4r,x5r,x6r,0.5*(x6+x6r)],[0,y1,y2,y3,y4,y5,y6,y1r,y2r,y3r,y4r,y5r,y6r,0.5*(y6+y6r)],[0,z1,z2,z3,z4,z5,z6,z1r,z2r,z3r,z4r,z5r,z6r,0.5*(z6+z6r)+5]
	ax.scatter(X,Y,Z,c = 'black',marker = 'o')
	
	line, = ax.plot([x1,x2], [y1,y2],[z1,z2], 'black', lw=1)
	line, = ax.plot([x2,x3], [y2,y3],[z2,z3], 'black', lw=1)
	line, = ax.plot([x3,x4], [y3,y4],[z3,z4], 'black', lw=2)
	line, = ax.plot([x4,x5], [y4,y5],[z4,z5], 'black', lw=1)
	line, = ax.plot([x5,x6], [y5,y6],[z5,z6], 'black'  , lw=1)

	line, = ax.plot([x1r,x2r], [y1r,y2r],[z1r,z2r], 'black', lw=1)
	line, = ax.plot([x2r,x3r], [y2r,y3r],[z2r,z3r], 'black', lw=1)
	line, = ax.plot([x3r,x4r], [y3r,y4r],[z3r,z4r], 'black', lw=2)
	line, = ax.plot([x4r,x5r], [y4r,y5r],[z4r,z5r], 'black', lw=1)
	line, = ax.plot([x5r,x6r], [y5r,y6r],[z5r,z6r], 'black'  , lw=1)

	line, = ax.plot([x6,x6r], [y6,y6r],[z6,z6r], 'black'  , lw=2)
	line, = ax.plot([x5,x5r], [y5,y5r],[z5,z5r], 'black'  , lw=1)
	
	print Tq[0]
	ax1.set_title("hip torque")
	ax1.set_xlim([0,10])
	ax1.set_ylim([-5,5])
	ax1.plot(t,Tq[0],'ro',markersize = 1)
	ax1.grid()

	ax2.set_title("knee torque")
	ax2.set_xlim([0,10])
	ax2.set_ylim([-10,10])	
	ax2.plot(t,Tq[1],'bo',markersize = 1)
	ax2.grid()

	ax3.set_title("hip angle")
	ax3.set_xlim([0,10])
	#ax3.set_ylim([-90,90])	
	ax3.plot(t,degrees(theta2),'bo',markersize = 1)
	ax3.grid()

	ax4.set_title("knee angle")
	ax4.set_xlim([0,10])
	#ax4.set_ylim([-90,90])	
	ax4.plot(t,degrees(theta3),'bo',markersize = 1)
	ax4.grid()

	ax.grid()
	plt.pause(0.000000000000000001)
	ax.cla()
