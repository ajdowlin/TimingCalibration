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

######## Configurables #########
inputFile = './PulsePulseCh23_1V_6kEvents.csv'

calibFilea_even = './Channel2_50158Events_EvenWindowCalib.csv'
calibFilea_odd = './Channel2_49840Events_OddWindowCalib.csv'
calibFileb_even = './Channel3_49913Events_EvenWindowCalib.csv'
calibFileb_odd = './Channel3_50085Events_OddWindowCalib.csv'

Cal = True #True if using calibrated x axis, False if uncalibrated
ClockSynced = False #True if not using relative measurement (use calibFilea)
TOA_low = 1000 #Omit events with too low of TOA
TOA_high = 10000 #Omit events with too high of TOA
numEvents = 5000 #Number of Events to use

################################

class Event:
    def __init__(self):
        self.windnum = []
        self.time = []
        self.data_cha = []
        self.data_chb = []
        
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
        EventsList[event_no].data_cha.append(fullCSV[row][4])
        EventsList[event_no].data_chb.append(fullCSV[row][5])
        row = row + 1
        if event_no%1000 == 0 and fullCSV[row][3] == 0:
            print("Reading Event ", event_no)
        
    return EventsList

def csvReader_sync(filename, numEvents):
    
    fullCSV = np.genfromtxt(filename, delimiter=',', skip_header=1)
    EventsList = [Event() for _ in range(numEvents)]
    event_no = 0
    row = 0                  
    while event_no < numEvents-1:
        event_no=int(fullCSV[row][1])
        EventsList[event_no].windnum.append(fullCSV[row][2])
        EventsList[event_no].time.append(fullCSV[row][3])
        EventsList[event_no].data_cha.append(fullCSV[row][4])
        row = row + 1
        if event_no%1000 == 0 and fullCSV[row][3] == 0:
            print("Reading Event ", event_no)
        
    return EventsList
    

def calcTOA(method, pulsex, pulsey, n,  peak = 500, perc = 0.7):

    #returns time  of threshold crossing using one of two methods:

    #inputs
    
    #method = 'fixed' uses a threshold of 300ADC counts
    #method = 'calc' uses a percentage of the calculated max for threshold
    #pulsex: uncalibrated x axis data of signal
    #pulsey: y axis data of signal
    #peak: set maximum used if method = 'fixed'
    #perc: percentage of peak used to find CFD point
    #n: number of points to use on each side of closest threshold crossing point for linear interp
    #ie, n=0 uses 2 points to do the linear interp, n=2 would use 6 points (2+2n)


    peakLoc = np.argmax(pulsey)

    if method == 'fixed':
        percPeak = peak * perc
    if method == 'calc':
        peakLeft = peakLoc - 2
        peakRight = peakLoc + 3
        calcpeak = np.average(pulsey[peakLeft:peakRight])
        percPeak = calcpeak * perc

    

    indexRight = np.searchsorted(pulsey[:peakLoc], percPeak) 
    indexLeft = indexRight - 1 
    x_linear = pulsex[(indexLeft-n):(indexRight+n+1)]
    y_linear = pulsey[(indexLeft-n):(indexRight+n+1)]

    slope, yintercept = np.polyfit(x_linear, y_linear, 1)
    TOA = (percPeak-yintercept)/slope

    return TOA


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


            
if ClockSynced == False:

    ch_adtEven = np.genfromtxt(calibFilea_even, delimiter=',', skip_header=1)
    ch_adtOdd = np.genfromtxt(calibFilea_odd, delimiter = ',', skip_header=1)
    ch_bdtEven = np.genfromtxt(calibFileb_even, delimiter=',', skip_header=1)
    ch_bdtOdd = np.genfromtxt(calibFileb_odd, delimiter = ',', skip_header=1)
    
    EventsList = csvReader(inputFile, numEvents)

    TOA_list = []
    
    for i in range(1, numEvents-1):
        startWindow = EventsList[i].windnum[0]
        x_old =  EventsList[i].time
        
        y_cha = EventsList[i].data_cha
        y_chb = EventsList[i].data_chb
        
        x_new_cha = fixTimeAxis(x_old, startWindow, ch_adtEven, ch_adtOdd)
        x_new_chb = fixTimeAxis(x_old, startWindow, ch_bdtEven, ch_bdtOdd)
        
        
        if Cal == True:
            
            TOA_cha = calcTOA('calc', x_new_cha, y_cha, 2,  peak = 400, perc = 0.71)
            TOA_chb = calcTOA('calc', x_new_chb, y_chb, 2,  peak = 400, perc = 0.71)
            
        if Cal == False:
                
            TOA_cha = calcTOA('calc', np.asarray(x_old)*100, y_cha, 2,  peak = 400, perc = 0.71)
            TOA_chb = calcTOA('calc', np.asarray(x_old)*100, y_chb, 2,  peak = 400, perc = 0.71)

        

        TOA = abs(TOA_cha - TOA_chb)
        if TOA > TOA_low and TOA < TOA_high:
            TOA_list.append(abs(TOA))
            


if ClockSynced == True:

    ch_adtEven = np.genfromtxt(calibFilea_even, delimiter=',', skip_header=1)
    ch_adtOdd = np.genfromtxt(calibFilea_odd, delimiter = ',', skip_header=1)

    EventsList = csvReader_sync(inputFile, numEvents)

    TOA_list = []
    
    for i in range(1, numEvents-1):
        startWindow = EventsList[i].windnum[0]
        x_old =  EventsList[i].time
        
        y_cha = EventsList[i].data_cha
        
        x_new_cha = fixTimeAxis(x_old, startWindow, ch_adtEven, ch_adtOdd)
        
        
        if Cal == True:
            
            TOA_cha = calcTOA('calc', x_new_cha, y_cha, 2,  peak = 400, perc = 0.71)
            
        if Cal == False:
                
            TOA_cha = calcTOA('calc', np.asarray(x_old)*100, y_cha, 2,  peak = 400, perc = 0.71)

        if abs(TOA_cha) > TOA_low and abs(TOA_cha) < TOA_high:
            TOA_list.append(abs(TOA_cha))
    

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

    

