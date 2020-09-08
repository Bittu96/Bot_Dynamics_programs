import sys
import time
import PyQt5 
from PyQt5 import *
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
import serial
import math
from math import *

def accdata(da):
    dat = str(da).strip()
    data = dat.replace('b','').replace('r','').replace('n','').replace("\\\\'",'').replace('\'','')
    x = ((180/pi) *(int(str(data[0:5]))  -20000)*0.1)
    y = ((180/pi) *(int(str(data[5:10])) -20000)*0.1)
    z = ((180/pi) *(int(str(data[10:15]))-20000)*0.1)
    return [x,y,z]


class joystickPlot(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(joystickPlot, self).__init__(parent)

        # Set up GUI configuration
        self.mainbox = QtGui.QWidget()
        self.setCentralWidget(self.mainbox)
        self.mainbox.setLayout(QtGui.QVBoxLayout())

        self.canvas = pg.GraphicsLayoutWidget()             # create GrpahicsLayoutWidget obejct  
        self.mainbox.layout().addWidget(self.canvas)

        self.label = QtGui.QLabel()                         # placeholder Qlabel object to display framerate
        self.mainbox.layout().addWidget(self.label)

     
        #  Set up plot
        self.analogPlot = self.canvas.addPlot(title='Real-time IMU plot')
        self.analogPlot.setYRange(-500,500)                # set axis range
        self.analogPlot.setXRange(-500,500)
        self.analogPlot.showGrid(x=True, y=True, alpha=0.5)
        x_axis = self.analogPlot.getAxis('bottom')
        y_axis = self.analogPlot.getAxis('left')
        
        x_axis.setLabel(text='x-axis reading')              # set axis labels
        y_axis.setLabel(text='y-axis reading')

        ticks = [[(0,'0'),(-500,'-500'),(500,'500')],[(-100,'-100'),(-200,'-200'),(-300,'-300'),(-400,'-400'),(100,'100'),(200,'200'),(300,'300'),(400,'400')]]
        x_axis.setTicks(ticks)                              # set ticks manually (first major-ticks than minor-ticks)
        y_axis.setTicks(ticks)
        
        self.drawplot = self.analogPlot.plot(pen='y')       # yellow line plot
        

        # initialize sensor data variables
        self.numPoints = 10                    # number of points that should be plottet at the same time stamp
        self.x = np.array([], dtype=int)        # create empty numpy array to store xAxis data
        self.y = np.array([], dtype=int)        # create empty numpy array to store yAxis data


        # initialize frame counter variables 
        self.counter = 0
        self.fps = 0.
        self.lastupdate = time.time()


        # set up image exporter (necessary to be able to export images)
        QtGui.QApplication.processEvents()
        self.exporter=pg.exporters.ImageExporter(self.canvas.scene())
        self.image_counter = 1
        
        
        # start updating
        self._update()



    def _framerate(self):
        now = time.time()   # current time since the epoch in seconds
        dt = (now-self.lastupdate)
        if dt <= 0:
            dt = 0.000000000001
        fps2 = 1.0 / dt
        self.lastupdate = now
        self.fps = self.fps * 0.9 + fps2 * 0.1
        tx = 'Mean Frame Rate:  {fps:.3f} FPS'.format(fps=self.fps )
        self.label.setText(tx)
        QtCore.QTimer.singleShot(1, self._update)
        self.counter += 1



    def _save_image(self):
        filename = 'img'+("%04d" % self.image_counter)+'.png'
        self.exporter.export(filename)
        self.image_counter += 1


    def _update(self):
        while (arduinoData.inWaiting()==0):                         #wait until there is data available
            pass # do nothing
        arduinoString = arduinoData.readline()  
        #print(arduinoString)                    #read the text line from serial port
        if sys.version_info >= (3, 0):
            arduinoString = arduinoString.decode('utf-8', 'backslashreplace')
            print(arduinoString) 
        dataArray = arduinoString.split(',')                        #split it into an array
        dataArray = accdata(arduinoString)
        
        xAxis = int(dataArray[0])                                   #convert first element to an integer
        yAxis = int(dataArray[1])                            #convert and flip yAxis signal
      
        self.x = np.append(self.x, xAxis)                           # append new data point to the array
        if self.x.size >= self.numPoints:                           # make sure that the size of the array includes the most recent data points
            self.x = np.append(self.x[1:self.numPoints],xAxis)
        else:
            self.x = np.append(self.x,xAxis) 

        self.y = np.append(self.y, yAxis)                           # append new data point to the array
        if self.y.size >= self.numPoints:                           # make sure that the size of the array includes the most recent data points
            self.y = np.append(self.y[1:self.numPoints],yAxis)
        else:
            self.y = np.append(self.y,yAxis)

        self.drawplot.setData(self.x, self.y)                       # draw current data set
        
        self._framerate()                                           # update framerate, see corresponding function

        # self._save_image()    # uncomment this to save each frame as an .png image in your current directory. Note that the framerate drops significantly by doing so

        

if __name__ == '__main__':
    arduinoData = serial.Serial("COM4", 57600);    # Creating our serial object named arduinoData, make sure that port and baud rate is set up according to your arduino data stream
    app = QtGui.QApplication(sys.argv)
    plot = joystickPlot()
    plot.show()
    sys.exit(app.exec_())