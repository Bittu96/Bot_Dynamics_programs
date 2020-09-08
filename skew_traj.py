import math
from math  import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt

fig = plt.figure()

T = 10
H = 5
X = 5
Y = 2

a = (8*H)/(T**3)
b = -(16*H)/(T**2)
c = (8*H)/(T)
d = 0

t = np.linspace(0,T,50)
tf = np.linspace(0,1,50)
Z_traverse = a*(t**3) + b*(t**2) + c*t
X_traverse = X*(((t/T) - (1/(2*pi))*cos((2*pi*t)/T - pi/2 )))
Y_traverse = Y*(((t/T) - (1/(2*pi))*cos((2*pi*t)/T - pi/2 )))
waist_traverse = abs(sin((pi/T)*t + pi/2))

heel_shift = 10*(1-cos(2*pi*tf/22))
heel_lift = 10*(sin(2*pi*tf/2))

r = 18
a = 1.5
b = 0
waist_traverse_c = sqrt(r**2 - (t-a)**2) + b

zmp = sign(sin((pi/0.5)*t))

plt.plot(t,Z_traverse)
plt.plot(t,X_traverse)
plt.plot(t,Y_traverse)
plt.plot(t,waist_traverse)
plt.plot(t,waist_traverse_c)
plt.plot(t,zmp,'black')
plt.plot(heel_shift,heel_lift,"red")

plt.xlim(0,20)
plt.ylim(0,20)

plt.grid()
plt.show()

'''
T = 5
x1 = np.linspace(0,T,50)
x2 = np.linspace(0,T/2,50)

p1 = 1+sin((2*pi*x1)/T - pi/2 )
p2 = 1+sin((2*pi*x2)/(T/2) - pi/2 )
p3 = (1+sin((2*pi*x2)/T - pi/2 )) + 0.5*(1+sin((2*pi*x2)/(T/2) - pi/2 ))

plt.plot(x1,p1)
plt.plot(x2,p2)
plt.plot(x2,p3)
'''