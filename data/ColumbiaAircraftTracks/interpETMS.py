import sys
import numpy as np

# from datetime import datetime
# refDate = datetime(1970,1,1,0,0,0)      #FACET refDate is refDate = datetime(1970,1,1,0,0,0)

import xml.etree.ElementTree as ET
tree = ET.parse('ETMS.20030201.134500-144500.Traf5')
root = tree.getroot()


# I confirm that the data in this ETMS file is decimal degrees.  Convert to HMS for input to FACET.
def convertDegreesToDegMinSec(val):
    Deg = int(np.floor(val))
    Min = int((val - Deg)*60.)
    Sec = int((val - Deg - Min/60.)*3600.)
    return [Deg, Min, Sec]


#Field10 appears to have flight plans?

### This comment block is to translate all of the aircraft in the air during the accident
# AircraftFromPaper = []
# ElementsFromPaper = root.iter('Flight')
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Just as a usage note, here's how to look through the different levels of the xml file
# flight.find('FlightPlans').find('FlightPlan').getchildren()
# where you notice that find gets only the first instance.  If you want subsequent instances, must use iter
# but beware because iter() will loop through all children to find the string


### This comment block is to capture just the elements in the list AircraftFromPaper
AircraftFromPaper = [
'COA282',
'CAA916',
'DAL1055',
'SWA333',
'COA1710',
'SKW3752',
'COA688',
'SKW3825',
'DAL2137']

ElementsFromPaper = []

for flight in root.iter('Flight'):
    AircraftID = flight.find('AircraftID').text
    if AircraftID in AircraftFromPaper:
        ElementsFromPaper.append(flight)

print ElementsFromPaper
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Detroit
# COA688 ->     FP_ROUTE IAH./.LFK.J29.ELD.J29.PXV..VHP..FWA.MIZAR3.DTW/1611

# Newark
# COA1710 ->     FP_ROUTE IAH./.LFK.J29.ELD.J29.MEM.J42.GVE.DYLIN2.EWR/1650

# Boston
# COA282 ->     FP_ROUTE IAH./.LFK.J29.ELD.J29.MEM.J42.RBV.J222.JFK.ORW3.BOS/1708

# Dallas
# CAA916 ->     FP_ROUTE HOU./.GIFFA.CQY5.DFW/1436
# SKW3825 ->     FP_ROUTE TYS./.ADMIT..SQS.CQY5.DFW/1443
# DAL1055 ->     FP_ROUTE PBI./.SRQ105014..AEX.CQY5.DFW/1444
# DAL2137 ->     FP_ROUTE MSY..WALDP..AEX.CQY5.DFW/0108
# SKW3752 ->     FP_ROUTE TLH./.CEW090069..CEW281044..CEW277049..AEX.CQY5.DFW/1437

# Nashville
# SWA333 ->     FP_ROUTE SAT./.WEMAR261027..LFK..LFK032027..EMG.J29.MEM.GHM4.BNA/1521

from Simulation.TJC import HaversineDistance as Haversine

from datetime import datetime
refDate = datetime(2003,2,1,0,0,0)      #The date of the Columbia accident

def secondsFromRefDate(curDateTime):
    dt = curDateTime - refDate    # Time delta since refDate
    return dt.total_seconds()

ACid = 0
aircraftRecord = {}
for flight in ElementsFromPaper:
    Callsign            = flight.find('AircraftID').text
    ACtype              = flight.find('AircraftType').text

    # if Callsign != 'SKW3825':
    #     continue

    fpVec               = []
    fpDatetimeVec       = []

    try:
        # curFlightPlan = flight.find('FlightPlans').find('FlightPlan')
        # curPlan = curFlightPlan.find('Field10').text

        # Get all flight plans
        fpVec           = [obj.text for obj in flight.find('FlightPlans').iter('Field10')]

        # Get all message times for flight plan changes (I don't know the difference between these two, but they are NOT same)
        fpTimeVec       = [obj.text for obj in flight.find('FlightPlans').iter('MessageTime')]
        fpDateVec       = [obj.text for obj in flight.find('FlightPlans').iter('MessageDate')]
        fpDatetimeVec   = [datetime.strptime(d + ' ' + t, '%Y-%m-%d %H:%M:%S') for (d,t) in zip(fpDateVec,fpTimeVec)]
    except:
        print 'NO FLIGHT PLAN -- SKIPPING'
        continue

    try:
        curTrackPoints  = flight.find('Track').find('Trackpoints')   # Only one Trackpoints, but has many Trackpoint children

        latVec          = [float(obj.text) for obj in curTrackPoints.iter('Latitude')]
        lonVec          = [float(obj.text) for obj in curTrackPoints.iter('Longitude')]

        filedSpeed      = [float(obj.text) for obj in curTrackPoints.iter('FiledSpeed')]

        timeVec         = [obj.text for obj in curTrackPoints.iter('Time')]
        dateVec         = [obj.text for obj in curTrackPoints.iter('Date')]

        # datetimeVec     = [datetime.strptime(d + ' ' + t, '%Y-%m-%d %H:%M:%S') for (d,t) in zip(dateVec,timeVec)]
        datetimeVec     = [secondsFromRefDate(datetime.strptime(d + ' ' + t, '%Y-%m-%d %H:%M:%S')) for (d,t) in zip(dateVec,timeVec)]
        flightLevelVec  = [int(obj.text) for obj in curTrackPoints.iter('Altitude')]

        # secondsFromMidnight = [secondsFromRefDate]
    except:
        print 'NO TRACKPOINTS -- SKIPPING'
        continue

    # if Callsign not in aircraftRecord:
    #     aircraftRecord[Callsign] = {'time':[], 'tracks':[], 'plans':[]}

    if Callsign in aircraftRecord:
        print "ERROR: {0} has multiple entries".format(Callsign)
        sys.exit()

    # aircraftRecord[Callsign] = {'time':datetimeVec, 'tracks':zip(latVec, lonVec), 'plans':fpVec}
    aircraftRecord[Callsign] = {'ACid':ACid, 'ACtype':ACtype, 'time':datetimeVec, \
                                'latVec':latVec, 'lonVec':lonVec, 'plans':fpVec, 
                                'filedSpeed':filedSpeed, 'flightLevelVec':flightLevelVec}
    ACid += 1


km_s2knots  = 1./0.000514444444


# Loss of signal is at about 13:59:32 UTC.  So let's say there's no danger before 13:55:00.  
# Similarly, let's say no danger after about two hours, 16:00:00.
# In fact, for the few planes we're considering here, they'll be out of the frame in under an hour
timeLo = secondsFromRefDate(datetime(2003,2,1,13,55,00))
timeHi = secondsFromRefDate(datetime(2003,2,1,15,00,00))

# Now interpolate on the aircraft record down to one-second intervals
secondsVec = np.linspace(timeLo, timeHi, timeHi - timeLo + 1) # Expand out all the seconds

trackTimeRecord = {}
# Preallocate 
for curTime in secondsVec:
    trackTimeRecord[int(curTime)] = []

for callSign in aircraftRecord:
    curAC = aircraftRecord[callSign]

    # Unpack the position info
    ACid = curAC['ACid']
    ACtype = curAC['ACtype']
    latVec = curAC['latVec']
    lonVec = curAC['lonVec']
    flightLevelVec = curAC['flightLevelVec']

    # Do the interpolations
    interpLat = np.interp(secondsVec, curAC['time'], latVec, left=0., right=0.)
    interpLon = np.interp(secondsVec, curAC['time'], lonVec, left=0., right=0.)
    interpFL = np.interp(secondsVec, curAC['time'], flightLevelVec, left=0., right=0.)

    # Calculate the ground speed
    ktasVec         = np.zeros_like(interpLon)     # Initialize to zero
    for ix in range(len(interpLat)-1):
        distKM = Haversine((interpLat[ix], interpLon[ix]), (interpLat[ix+1], interpLon[ix+1]))
        ktasVec[ix] = (distKM / 1.) * km_s2knots              # delta_t is one second, so km/s to knots
    ktasVec[-1] = ktasVec[-2] # Set the last value to be equal to the second-to-last

    # Store them by times
    for tx in range(len(secondsVec)):
        curTime = int(secondsVec[tx])
        latDec = interpLat[tx]
        lonDec = interpLon[tx]
        curLine = "{0} {1} {2:.6f} {3:.6f} {4:.2f} {5:.2f}".format(ACid, ACtype, latDec, lonDec, int(interpFL[tx]), ktasVec[tx])
        if latDec > 0:
            trackTimeRecord[curTime].append(curLine)

# So here's the filed speed.  Multiple elements, not just [0]
# flight.find('FlightPlans').getchildren()[0].find('FiledSpeed').text

# Now dump to file
# Open the file and get the dict keys
outFileName = 'HighRiskETMSInterp.txt'
try:
    outFile = open(outFileName,'w')
except:
    print 'file failed to open'

sortedKeys = sorted(trackTimeRecord.keys())
for curTime in sortedKeys:
    curDict = trackTimeRecord[curTime]

    outFile.write("{0}\n".format(curTime))

    for line in curDict:
        # print curDict
        outFile.write(line + '\n')
    outFile.write('\n')
outFile.close()


# Dump the callsign keys
outFileName = 'HighRiskETMSInterp_CallsignKeys.txt'
try:
    outFile = open(outFileName,'w')
except:
    print 'file failed to open'

for callSign in aircraftRecord:
    curDict = aircraftRecord[callSign]

    outFile.write("{0} {1} {2}\n".format(curDict['ACid'], callSign, curDict['ACtype']))

    # for line in curDict:
    #     # print curDict
    #     outFile.write(line + '\n')
    # outFile.write('\n')
outFile.close()












