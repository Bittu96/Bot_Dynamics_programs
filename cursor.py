import pyautogui, sys
import serial
from pynput.mouse import *

mouse = Controller()
baudrate = 9600
ser=serial.Serial('COM13',baudrate)
print(ser.name)
ser.baudrate = baudrate

def accdata(data):
    dataX = str(data.decode('utf-8', 'backslashreplace'))[0:8]
    dataY = str(data.decode('utf-8', 'backslashreplace'))[9:17]
    #print(dataY)
    global valx,valy
    valx= valy= 0
    for i in range(0,8):
        valx += (ord(list(dataX)[i])-48)* 10**(4-i)
        valy += (ord(list(dataY)[i])-48)* 10**(4-i)
    return valx-20000,valy-20000

k = 0
while True:
    #print((ser.readline()).decode('utf-8', 'backslashreplace'))
    acc = accdata(ser.readline())
    #print(int(acc[0])+500, int(acc[1])+500)
    #print(mouse.position)
    mouse.position = (-int(acc[1])*10 +500,int(acc[0])*10 +500)
    #print(ser.readline())
    # if k&10 == 0:
    #     print('bbhehe')
    #     pyautogui.moveTo(int(acc[1])*10 +500, int(acc[0])*10 +500)
    # k+=1
