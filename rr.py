import math
from math import *
import numpy
from numpy import *
import time
import pygame
from pygame import *
import sys, pygame, math, os, random
from pygame.locals import *

'''
import sympy
from sympy import *

theta1 = Symbol("theta1")
theta2 = Symbol("theta2")

theta1d = Symbol("theta1d")
theta2d = Symbol("theta2d")

theta1dd = Symbol("theta1dd")
theta2dd = Symbol("theta2dd")
'''
#\\\\\\\\\\\\\  values\\\\\\\\\\\\\\\\\\\
pygame.init()
size=width,height=800,600
screen=pygame.display.set_mode(size)
pygame.display.set_caption("RR bot Dynamics Simulation")

SimulationRunning = True

theta1i = pi/2 + pi/3 #radians
theta2i = 0 + pi/3

theta1di = 0 #radians/sec
theta2di = 0

theta1ddi = 0 #radians/sec2
theta2ddi = 0

#\\\\\\\\\\\\\  PARAMETERS \\\\\\\\\\\\\\\\
m1 = 5      #kg
m2 = 5      #kg
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

height_step = 0.1 #m
T_step = 3.0        #sec
time_interval = 0.01
theta1_track = height_step*sin(2*pi*time_interval*(1/T_step))
theta2_track = height_step*sin(2*pi*time_interval*(1/T_step))

####################################### VISUAL #############################################################################
####################################### VISUAL #############################################################################
####################################### VISUAL #############################################################################

time = 0
scalingFactor = sf = 500

while SimulationRunning:

  bgcolour=(0,0,0) 
  screen.fill(bgcolour)
  
  theta1_track = height_step*sin(2*pi*time*(1/T_step))
  theta2_track = height_step*sin(2*pi*time*(1/T_step))
  
  theta1 = theta1i + theta1_track
  theta2 = theta2i - 2*theta2_track
  
  theta1d = theta1di + (2*pi*(1/T_step))*height_step*cos(2*pi*time*(1/T_step))
  theta2d = theta2di + (2*pi*(1/T_step))*height_step*cos(2*pi*time*(1/T_step))

  theta1dd = theta1ddi - height_step*sin(2*pi*time*(1/T_step))
  theta2dd = theta2ddi - height_step*sin(2*pi*time*(1/T_step))

  torqueCalc(theta1,theta2,theta1d,theta2d,theta1dd,theta2dd)

  joint1 = [400,300]
  joint2 = [ int( joint1[0] + l1*sf*sin(theta1) )               , int( joint1[1] + l1*sf*cos(theta1) ) ]
  endeff = [ int( joint2[0] + l2*sf*sin(theta2+theta1) )        , int( joint2[1] + l2*sf*cos(theta2+theta1) ) ]

  pygame.draw.line(screen, (255,0,0), joint1, joint2, 2)
  pygame.draw.line(screen, (255,0,0), joint2, endeff, 2)
 
  for event in pygame.event.get(): 
    if event.type==pygame.QUIT:sys.exit()
  pygame.draw.circle(screen, (255,0,0)  , joint1, 10,2)
  pygame.draw.circle(screen, (255,0,0)  , joint2, 10,2)
  pygame.draw.circle(screen, (255,0,0)  , endeff, 10,2) 
  pygame.display.flip();pygame.time.delay(10)
  pygame.display.update()
  time += 0.01
