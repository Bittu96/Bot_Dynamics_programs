import math
from math import *
import numpy
from numpy import *
import time
import pygame
from pygame import *
import sys, pygame, math, os, random
from pygame.locals import *

pygame.init()
size=width,height=1800,500
screen=pygame.display.set_mode(size)
pygame.display.set_caption("inverted_pendulum_balancing")

th = pi/10
theta = pi - th

g = 9.81
m = 5
l = 300
t = 0
base = [100,370]
acc = 0
the = 0
pos=[(base[0],base[1]), ( base[0] + l*sin(theta),(base[1] + l*cos(theta)))]  

training_data = [[],[],[]]
t = 0
av = 0
while t<10:
  #
  #acc = 5*g*tan(th)
  th += ( (0.5/l) * (g*sin(th) - acc*cos(th)) )
  training_data[1].append(th/pi)
  training_data[2].append(av-th)
  av=th 
  key=pygame.key.get_pressed()
  if key[pygame.K_UP]:   
    acc+=0.8
    training_data[0].append(1)
  elif key[pygame.K_DOWN]: 
    acc-=0.8
    training_data[0].append(-1)
  else:
    acc=(acc*0.99)
    training_data[0].append(0)
  base[0] += int(acc)
  print mat(training_data).T#degrees(th),acc
  theta = pi - th
  pos=[(base[0],base[1]), (int(base[0] + l*sin(theta)),int(base[1] + l*cos(theta)))]   
  the -= 0.015*acc
  t += 0.01
  #
  bgcolour=(0,0,0) 
  screen.fill(bgcolour)
  pygame.draw.line(screen, (255,0,0), (0,400), (1600,400), 2)
  pygame.draw.line(screen, (255,0,0), pos[0], pos[1], 2)
  pygame.draw.line(screen, (0,255,255), (pos[0][0] - 30*sin(the),pos[0][1] - 30*cos(the)), (pos[0][0] + 30*sin(the),pos[0][1] + 30*cos(the)), 2)
  pygame.draw.line(screen, (0,255,255), (pos[0][0] - 30*sin(the+pi/1.5),pos[0][1] - 30*cos(the + pi/1.5)), (pos[0][0] + 30*sin(the+pi/1.5),pos[0][1] + 30*cos(the+pi/1.5)), 2)
  pygame.draw.line(screen, (0,255,255), (pos[0][0] - 30*sin(the-pi/1.5),pos[0][1] - 30*cos(the - pi/1.5)), (pos[0][0] + 30*sin(the-pi/1.5),pos[0][1] + 30*cos(the-pi/1.5)), 2)
  
  for event in pygame.event.get(): 
    if event.type==pygame.QUIT:sys.exit()
  pygame.draw.circle(screen, (255,0,0), pos[1], 30,5) 
  pygame.draw.circle(screen, (0,255,255), pos[0], 30,3)
  pygame.display.flip();pygame.time.delay(10)
  pygame.display.update()

#print mat(training_data[0]).T
#print mat(training_data[1]).T
print "TRAINING...."

l = 1

def sigm(x, deriv=False):
  if(deriv==True):
    return (x.T*(1-x))
  return 1/(1+ exp(-x))

w0 = (1,20)
w1 = (20,1)

wt0 = numpy.random.uniform(low=0, high=1, size=w0)
wt1 = numpy.random.uniform(low=0, high=1, size=w1)

x = array(training_data[1]).T
for i in range(0,50):
  print i
  x = array(training_data[1]).T[i]
  y = array(training_data[0]).T[i]
 
  for j in xrange(20000):
    l0 = mat(x)
    l1 = sigm(dot(l0, wt0))
    l2 = sigm(dot(l1, wt1))

    l2_error = y - l2
    if (j % 10000) == 0:
      print 'Error:' + str(mean(abs(l2_error)))

    l2_delta = l2_error*sigm(l2, deriv=True)
    l1_error = l2_delta.dot(wt1.T)
    l1_delta = l1_error*sigm(l1, deriv=True)

    wt1 += l*l1.T.dot(l2_delta)
    wt0 += l*l0.T.dot(l1_delta)

print 'After training'  
print l2
print wt0
print wt1
'''
pygame.init()
size=width,height=1800,1600
screen=pygame.display.set_mode(size)
pygame.display.set_caption("inverted_pendulum_balancing")

while 1:
  th += ( (0.5/l) * (g*sin(th) - acc*cos(th)) )
  
  l0 = mat(th)
  l1 = sigm(dot(l0, wt0))
  l2 = sigm(dot(l1, wt1))
  
  if   l2 > 0:   
    acc+=0.8
  elif l2 < 0: 
    acc-=0.8
  else:
    acc=(acc*0.99)
  
  base[0] += int(acc)
  print acc
  theta = pi - th
  pos=[(base[0],base[1]), (int(base[0] + l*sin(theta)),int(base[1] + l*cos(theta)))]   
  the -= 0.015*acc
  t += 0.01
  #
  bgcolour=(0,0,0) 
  screen.fill(bgcolour)
  pygame.draw.line(screen, (255,0,0), (0,400), (800,400), 2)
  pygame.draw.line(screen, (255,0,0), pos[0], pos[1], 2)
  pygame.draw.line(screen, (0,255,255), (pos[0][0] - 30*sin(the),pos[0][1] - 30*cos(the)), (pos[0][0] + 30*sin(the),pos[0][1] + 30*cos(the)), 2)
  pygame.draw.line(screen, (0,255,255), (pos[0][0] - 30*sin(the+pi/1.5),pos[0][1] - 30*cos(the + pi/1.5)), (pos[0][0] + 30*sin(the+pi/1.5),pos[0][1] + 30*cos(the+pi/1.5)), 2)
  pygame.draw.line(screen, (0,255,255), (pos[0][0] - 30*sin(the-pi/1.5),pos[0][1] - 30*cos(the - pi/1.5)), (pos[0][0] + 30*sin(the-pi/1.5),pos[0][1] + 30*cos(the-pi/1.5)), 2)
  
  for event in pygame.event.get(): 
    if event.type==pygame.QUIT:sys.exit()
  pygame.draw.circle(screen, (255,0,0), pos[1], 30,5) 
  pygame.draw.circle(screen, (0,255,255), pos[0], 30,3)
  pygame.display.flip();pygame.time.delay(10)
  pygame.display.update()
'''