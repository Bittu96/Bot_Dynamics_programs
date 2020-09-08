import serial
from math import *
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure()
ax = Axes3D(fig)

baudrate = 57600
ser=serial.Serial('COM13',baudrate)
print(ser.name)
ser.baudrate = baudrate


def accdata(data):
    data = str(data.decode('utf-8', 'backslashreplace'))
    print(data)
    x = int(data[0:5])-20000
    y = int(data[6:11])-20000
    z = int(data[12:17])-20000
    return [x,y,z]


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


def dh(t,d,a,al):
    fin = T(0,0,t,0,0,0)*T(0,0,0,0,0,d)*T(0,0,0,a,0,0)*T(al,0,0,0,0,0)
    return fin



while True:
    #print((ser.readline()).decode('utf-8', 'backslashreplace'))
    #print(accdata(ser.readline()))
    print((ser.readline()))
    #print(acc)
