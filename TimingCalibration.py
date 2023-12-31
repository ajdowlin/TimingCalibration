import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import mplhep

######## Configurables #########

inputFile = '20kEvents_1V_1ns_Gauss_ch0.csv'

plotHist = True #show histogram of threshold crossings
plotSpread = True #show distribution of measured delays
Exaggerate = True #True to 'exaggerate' measured time intervals
exaggerate_percent = 0.1 #how much to 'exaggerate' measured time intervals if Exaggerate = True
output = True #if true, creates csv output file
channel = 1 #specify channel used

################################


class Event:
    def __init__(self):
        self.windnum = []
        self.time = []
        self.data_ch0 = []

        
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
        row = row + 1
        if event_no%10000 == 0 and fullCSV[row][3] == 0:
            print("Reading Event ", event_no)
        
    return EventsList, numEvents

def calcTOA(method, pulsex, pulsey, n,  peak = 800, perc = 0.7):

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
        peakLeft = peakLoc - 3
        peakRight = peakLoc + 3
        peak = np.average(pulsey[peakLeft:peakRight])
        percPeak = peak * perc

    

    indexRight = np.searchsorted(pulsey[:peakLoc], percPeak)  
    indexLeft = indexRight - 1 
    x_linear = pulsex[(indexLeft-n):(indexRight+n)]
    y_linear = pulsey[(indexLeft-n):(indexRight+n)]

    slope, yintercept = np.polyfit(x_linear, y_linear, 1)
    TOA = (percPeak-yintercept)/slope

    return TOA

################################
#Load Events into data structure
################################

print("Loading events...")
EventsList, numEvents = csvReader(inputFile, 20000)
events_used = range(1, numEvents-1) #omit first and last event

#####################################################
#Find Sample Crossings for each event, fill histogram
#####################################################

TOA = []
dtOdd = np.zeros(64, dtype = int)
dtEven = np.zeros(64, dtype = int)

for i in events_used:
    if i%10000 == 0:
        print("Finding threshold cossing for event ", i)
    x, y, windnum = EventsList[i].time, EventsList[i].data_ch0, EventsList[i].windnum
    TOA = calcTOA('fixed', x, y, 1)

    if TOA > len(windnum): ## omit some code-breaking outliers
        continue
    
    win = windnum[int(TOA)]
    sample_num = TOA%64
    
#account for even and odd window differences

    if win%2 == 0:
        dtEven[int(sample_num)] = dtEven[int(sample_num)] + 1

    else:
        dtOdd[int(sample_num)] = dtOdd[int(sample_num)] + 1

lenOdd = np.sum(dtOdd)
lenEven = np.sum(dtEven)

print("Number of samples located in Odd windows: ", lenOdd)
print("Number of samples located in Even windows: ", lenEven)

#Scale histogram to known length of window


if Exaggerate == True:
    min_crossings_even = np.min(dtEven)
    min_crossings_odd = np.min(dtOdd)

    toRemove_even = np.min(dtEven)*exaggerate_percent
    toRemove_odd = np.min(dtOdd)*exaggerate_percent

    dtEven = dtEven - toRemove_even
    dtOdd = dtOdd - toRemove_odd

dtOdd_scale = (dtOdd/np.sum(dtOdd)) * 6400
dtEven_scale = (dtEven/np.sum(dtEven)) * 6400

dtOdd_DF = pd.DataFrame(dtOdd_scale)
dtEven_DF = pd.DataFrame(dtEven_scale)

dt = list(dtEven)+list(dtOdd)

if plotHist == True:
    mplhep.style.use("LHCb2")
    fig, axes = plt.subplots()
    axes.set_xlabel("Threshold crossing location", loc = 'center', fontsize = 15)
    axes.bar(range(128), dt)
    fig.set_size_inches(10,7)
    plt.show()
        
if plotSpread == True:
    mplhep.style.use("LHCb2")
    fig, axes = plt.subplots()
    axes.set_xlabel("True delay", loc = 'center', fontsize = 15)
    axes.hist(dt, 15)
    fig.set_size_inches(10,7)
    plt.show()

if output == True:
    outString_even = "Ch" + str(channel) + "_" + str(lenEven) + "Events_EvenWinCalib"
    outString_odd = "Ch" + str(channel) + "_" + str(lenOdd) + "Events_OddWinCalib"
    if Exaggerate  == True:
        outString_even = outString_even + "_" + str(exaggerate_percent) + "exaggerate"
        outString_odd = outString_odd + "_" + str(exaggerate_percent) + "exaggerate"
        
        
    dtEven_DF.to_csv(outString_even + ".csv")
    dtOdd_DF.to_csv(outString_odd + ".csv")
    

