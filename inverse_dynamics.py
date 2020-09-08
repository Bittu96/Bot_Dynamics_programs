from math import sin,cos
from numpy import array,mat,ravel
import sys 
from termcolor import colored, cprint 

####################---------Inverse Dynamics of RRR
# text = colored('Hello, World!', 'red', attrs=['reverse', 'blink']) 
# print(text)
# print(colored('Hello, World!', 'green', 'on_red')) 

L1=0.18
L2=0.18
L3=0.05
M1=1.5
M2=1.5
M3=1.0
spec=[L1, L2, L3, M1, M2, M3]

def TorqueCalc(x):
    global spec
    theta1dd,theta2dd,theta3dd = 3,-11,-1
    qdd = mat([ [theta1dd], [theta2dd], [theta3dd] ])

    L1 = spec[0]
    L2 = spec[1]
    L3 = spec[2]
    M1 = spec[3]
    M2 = spec[4]
    M3 = spec[5]
    g = 9.8

    #print('hehe :',x[0])

    ############################-------INERTIAL TERM---------#################

    I11 = ( 
        (M1+M2+M3)*L1**2 +
        (M2+M3)*L2**2 +
        M3*L3**2 + 
        2*(M2+M3)*L1*L2*cos(x[5]) + 
        2*M3*L2*L3*cos(x[6]) + 
        2*M3*L3*L1*cos(x[5]+ x[6])
        )
    I12 = (
        (M2+M3)*L2**2 + 
        (M2+M3)*L1*L2*cos(x[5]) + 
        M3*L3**2 + 
        2*M3*L2*L3*cos(x[6]) + 
        M3*L3*L1*cos(x[5]+ x[6])
        )
    I13 = (
        M3*L3**2 + 
        M3*L2*L3*cos(x[6]) + 
        M3*L3*L1*cos(x[5]+ x[6])
        )
    I21 = (
        (M2+M3)*L2**2 + 
        (M2+M3)*L1*L2*cos(x[5]) + 
        M3*L3**2 + 
        2*M3*L2*L3*cos(x[6]) + 
        M3*L3*L1*cos(x[5]+ x[6])
        )
    I22 = ( 
        (M2+M3)*L2**2 + 
        M3*L3**2 + 
        2*M2*L2*L3*cos(x[6])
        )
    I23 = ( 
        M3*L3**2 + 
        M3*L2*L3*cos(x[6])
        )
    I31 = (
        M3*L3**2 + 
        M3*L2*L3*cos(x[6]) + 
        M3*L3*L1*cos(x[5]+ x[6])
        )
    I32 = M3*L3**2 + M3*L2*L3*cos(x[6])
    I33 = M3*L3**2
    
    I  = mat([ [I11, I12, I13], 
               [I21, I22, I23], 
               [I31, I32, I33] ])

    ############################-------INERTIAL TERM---------#################
    ############################-------CORIOLIS TERM---------#################
    
    c1 = (
        -(2*(M2+M3)*L1*L2*sin(x[5]) + 2*M3*L3*L1*sin(x[5] + x[6]))*x[7]*x[8] - 
        (2*M3*L2*L3*sin(x[6]) + 2*M3*L3*L1*sin(x[5] + x[6]))*x[8]*x[9] - 
        (2*M3*L2*L3*sin(x[6]) + 2*M3*L3*L1*sin(x[5] + x[6]))*x[9]*x[7] - 
        ((M2+M3)*L1*L2*sin(x[5]) + M3*L3*L1*sin(x[5] + x[6]))*(x[8])**2 - 
        (M3*L2*L3*sin(x[6]) + M3*L3*L1*sin(x[5] + x[6]))*(x[9])**2
        )
    c2 = (
        -((M2+M3)*L1*L2*sin(x[5]) + 
        M3*L3*L1*sin(x[5] + 
        x[6])*x[7]*x[8]) - 
        (2*M3*L2*L3*sin(x[6])*x[8]*x[9]) - 
        (2*M3*L2*L3*sin(x[6]) + 
        M3*L3*L1*sin(x[5] + x[6]))*x[7]*x[9] - 
        (M3*L2*L3*sin(x[6])*(x[9])**2)
        )
    c3 = (
        -M3*L2*L3*sin(x[6])*x[8]*x[9] - 
        (M3*L2*L3*sin(x[6]) + 
        M3*L3*L1*sin(x[5] + 
        x[6]))*x[7]*x[8] - 
        M3*L3*L1*sin(x[5] + x[6])*x[7]*x[8]
        )
    
    C = mat([ [c1], 
              [c2], 
              [c3] ])

    ############################-------CORIOLIS TERM---------#################
    #############################-------GRAVITY TERM---------#################
    
    g1 = (
        -(M1+M2+M3)*g*L1*sin(x[4]) - 
        (M2+M3)*g*L2*sin(x[4] + x[5]) - 
        M3*g*L3*sin(x[4] + x[5] + x[6])
        )
    g2 = (
        -(M2+M3)*g*L2*sin(x[4] + x[5]) - 
        M3*g*L3*sin(x[4] + x[5] + x[6])
        )
    g3 = -M3*g*L3*sin(x[4] + x[5] + x[6])
    
    G = mat([ [g1], 
              [g2], 
              [g3] ])
    
    #############################-------GRAVITY TERM---------#################

    Torques = I*qdd + C + G 
    return Torques


x = array([0,1,0,0,0,0,0,0,0,0])
values = TorqueCalc(x)
print('Torques :',ravel(values))

th1,th2,th3 = 0,0,0
x1=L1*sin(th1)
y1=L1*cos(th1)
x2=L1*sin(th1)+L2*sin(th1+th2)
y2=L1*cos(th1)+L2*cos(th1+th2)
x3=L1*sin(th1)+L2*sin(th1+th2)+L3*sin(th1+th2+th3)
y3=L1*cos(th1)+L2*cos(th1+th2)+L3*cos(th1+th2+th3)