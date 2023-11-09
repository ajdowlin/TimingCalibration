import csv
import numpy as np
import math
import scipy
from iminuit import Minuit, cost
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import mplhep
import os

class Event:
    def __init__(self):
        self.windnum = []
        self.time = []
        self.data_ch0 = []
        self.data_ch3 = []
        
def csvReader(filename, numEvents):

    #Create datastructure from csv file
    
    fullCSV = np.genfromtxt(filename, delimiter=',', skip_header=1)
    EventsList = [Event() for _ in range(numEvents)]
    event_no = 0
    row = 0                  
    while event_no < numEvents-1:
        event_no=int(fullCSV[row][1])
        EventsList[event_no].windnum.append(fullCSV[row][2])
        EventsList[event_no].time.append(fullCSV[row][3])
        EventsList[event_no].data_ch0.append(fullCSV[row][4])
        EventsList[event_no].data_ch3.append(fullCSV[row][5])
        row = row + 1
        if event_no%1000 == 0 and fullCSV[row][3] == 0:
            print("Reading Event ", event_no)
        
    return EventsList, numEvents

def calcTOA(method, pulsex, pulsey, n,  peak = 500, perc = 0.7):

    #returns location of threshold crossing using one of two methods:
    #method = 'fixed' uses a threshold of 300ADC counts
    #method = 'calc' uses a percentage of the calculated max for threshold

    peakLoc = np.argmax(pulsey)

    if method == 'fixed':
        percPeak = peak * perc
    if method == 'calc':
        peakLeft = peakLoc - 2
        peakRight = peakLoc + 3
        peak = np.average(pulsey[peakLeft:peakRight])
        percPeak = peak * perc

    

    indexRight = np.searchsorted(pulsey[:peakLoc], percPeak) 
    indexLeft = indexRight - 1 
    x_linear = pulsex[(indexLeft-n):(indexRight+n+1)]
    y_linear = pulsey[(indexLeft-n):(indexRight+n+1)]

    slope, yintercept = np.polyfit(x_linear, y_linear, 1)
    TOA = (percPeak-yintercept)/slope

    return TOA


def trimGauss(x, y, xl, xr):

    #Trim event so that only pulse is used
    
    maxLoc = np.argmax(y)
    x_trim = x[(maxLoc-xl):(maxLoc+xr)]
    y_trim = y[(maxLoc-xl):(maxLoc+xr)]

    return x_trim, y_trim

def fixTimeAxis(x, startWin, dtEven, dtOdd):
    totalWindows = int(len(x)/64)
    dtEven_time = dtEven[:, 1]
    dtOdd_time = dtOdd[:, 1]

    x_new = [0]
    
    if startWin%2 == 0:
        dt = np.asarray(list(dtEven_time) + list(dtOdd_time))
        for i in range(int(totalWindows/2)):
            for val in dt:
                x_last = x_new[len(x_new)-1]
                x_new.append(val+x_last)
            
            
            

    if startWin%2 == 1:
        dt = np.asarray(list(dtOdd_time) + list(dtEven_time))
        for i in range(int(totalWindows/2)):
            for val in dt:
                x_last = x_new[len(x_new)-1]
                x_new.append(val+x_last)
       
                 

                
    return x_new
            
            
        
#EventsList, numEvents = csvReader('2000Events_newCables.csv', 500)
EventsList, numEvents = csvReader('PulsePulse_Gauss_1kEvents.csv', 1000)
ch_0dtEven = np.genfromtxt('Channel0_EvenWindowCalibration', delimiter=',', skip_header=1)
ch_0dtOdd = np.genfromtxt('Channel0_OddWindowCalibration', delimiter = ',', skip_header=1)
ch_3dtEven = np.genfromtxt('Channel3_EvenWindowCalibration', delimiter=',', skip_header=1)
ch_3dtOdd = np.genfromtxt('Channel3_OddWindowCalibration', delimiter = ',', skip_header=1)


TOA_list = []
Cal = False

for i in range(1, numEvents-1):
    startWindow = EventsList[i].windnum[0]
    x_old =  EventsList[i].time
    
    y_ch0 = EventsList[i].data_ch0
    y_ch3 = EventsList[i].data_ch3
    
    x_new_ch0 = fixTimeAxis(x_old, startWindow, ch_0dtEven, ch_0dtOdd)
    x_new_ch3 = fixTimeAxis(x_old, startWindow, ch_3dtEven, ch_3dtOdd)

    #x_ch0, y_ch0 = trimGauss(x_new_ch0, y_ch0, 30, 30)
    #x_ch3, y_ch3 = trimGauss(x_new_ch3, y_ch3, 30, 30)

    if Cal == True:
    
        TOA_ch0 = calcTOA('calc', x_new_ch0, y_ch0, 2,  peak = 400, perc = 0.71)
        TOA_ch3 = calcTOA('calc', x_new_ch0, y_ch3, 2,  peak = 400, perc = 0.71)
        
    if Cal == False:
        
        TOA_ch0 = calcTOA('calc', np.asarray(x_old)*100, y_ch0, 2,  peak = 400, perc = 0.71)
        TOA_ch3 = calcTOA('calc', np.asarray(x_old)*100, y_ch3, 2,  peak = 400, perc = 0.71)

        

    TOA = (TOA_ch0 - TOA_ch3)
    TOA_list.append(TOA)
    #if TOA > 10 and TOA < 13:
     #   TOA_list.append(TOA)

print("TOA Mean: ", np.mean(TOA_list))
print("TOA STD: ", np.std(TOA_list))
#print(x_new_ch0)

mplhep.style.use("LHCb2")
fig, axes = plt.subplots()
axes.set_xlabel("ps", loc = 'center', fontsize = 20)
#axes.plot(EventsList[62].time, EventsList[62].data_ch3)
#axes.plot(EventsList[62].time, EventsList[62].data_ch0)
axes.hist(np.asarray(TOA_list), 50)
fig.set_size_inches(10,7)
plt.show()

    

