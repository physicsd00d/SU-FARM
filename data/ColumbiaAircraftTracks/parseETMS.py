import sys
import numpy as np

from datetime import datetime
refDate = datetime(1970,1,1,0,0,0)      #FACET refDate is refDate = datetime(1970,1,1,0,0,0)

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


trackTimeRecord = dict()
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

        # # Format for lat lon are to be six-digit integers, lon is west (so positive)
        # latVec = [str(obj*1e4).split('.')[0] for obj in latVec]
        # lonVec = [str(abs(obj)*1e4).split('.')[0] for obj in lonVec]

        # Format for lat lon are to be six-digit integers, lon is west (so positive)
        latVec          = [convertDegreesToDegMinSec(obj) for obj in latVec] # This is list of 3-element lists, floats
        lonVec          = [convertDegreesToDegMinSec(abs(obj)) for obj in lonVec]

        timeVec         = [obj.text for obj in curTrackPoints.iter('Time')]
        dateVec         = [obj.text for obj in curTrackPoints.iter('Date')]

        datetimeVec     = [datetime.strptime(d + ' ' + t, '%Y-%m-%d %H:%M:%S') for (d,t) in zip(dateVec,timeVec)]
        flightLevelVec  = [obj.text for obj in curTrackPoints.iter('Altitude')]

        # # Not sure how to use airspeed because FACET like groundspeed
        # airspeedVec     = [obj.text for obj in curTrackPoints.iter('IndicatedAirspeed')]

    except:
        print 'NO TRACKPOINTS -- SKIPPING'
        continue

    # These are not given
    heading = '0'             # I think this doesn't matter
    curCenter = ' '         # We'll see if we need this one
    curSector = ' '         # NOPE!  You can totally get away without these :)

    # curPlanIX = -1  # Initialize this
    for (date, lat, lon, FL) in zip(datetimeVec, latVec, lonVec, flightLevelVec):

        # Which flight plan to use?
        curPlan = ''
        truthList = [(date <= obj) for obj in fpDatetimeVec]
        if (True in truthList):
            curPlanIX = [(date <= obj) for obj in fpDatetimeVec].index(True) - 1
        else:
            # Use last index
            curPlanIX = len(fpDatetimeVec)-1

        if (curPlanIX) == -1:
            # Use the first plan because we haven't hit the first update yet
            curPlan = fpVec[0]
        else:
            # The first time this gets hit will still be using the 0th flight plan but we'll be in the updating regime
            curPlan = fpVec[curPlanIX]


        latStr = '{0:02}{1:02}{2:02}'.format(lat[0], lat[1], lat[2])
        lonStr = '{0:02}{1:02}{2:02}'.format(lon[0], lon[1], lon[2])

        infoString1 = 'TRACK ' + Callsign + ' ' + ACtype + ' ' + latStr + ' ' + \
                lonStr + ' ' + '1' + ' ' + str(FL) + ' ' + heading + ' ' + curCenter + \
                ' ' + curSector
        infoString2 = '    FP_ROUTE ' + curPlan

        curDict = {'CallSign': Callsign,
           'Aircraft': ACtype,
           'infoString2': infoString2,
           'infoString1': infoString1}

        # Add the key to the dictionary or append to existing key
        if date not in trackTimeRecord:
            trackTimeRecord[date] = []
            trackTimeRecord[date].append(curDict)
        else:
            trackTimeRecord[date].append(curDict)



# fptimevec = [obj.find('CoordinationFixTime').text for obj in flight.find('FlightPlans').iter('FlightPlan')] #Departure Time

# [obj.text for obj in flight.find('FlightPlans').iter('FiledSpeed')]
# [obj.text for obj in flight.find('FlightPlans').iter('AssignedAltitude')]
# [obj.text for obj in curTrackPoints.iter('Altitude')]
# [obj.text for obj in curTrackPoints.iter('IndicatedAirspeed')]
# [obj.text for obj in curTrackPoints.iter('TrueAirspeed')]   #zeros
# [obj.text for obj in curTrackPoints.iter('Course')]   #zeros


# Open the file and get the dict keys
outFileName = 'TRX_Columbia_HighRisk003'
try:
    outFile = open(outFileName,'w')
except:
    print 'file failed to open'

sortedKeys = sorted(trackTimeRecord.keys())
for curTime in sortedKeys:
    curDict = trackTimeRecord[curTime]
    TRACK_TIME = int((curTime-refDate).total_seconds())
    outFile.write('TRACK_TIME ' + str(TRACK_TIME) + '\n')

    for flight in curDict:
        # print curDict
        outFile.write(flight['infoString1'] + '\n')
        outFile.write(flight['infoString2'] + '\n\n')
outFile.close()


# Make explicitly sure track lat lon in degrees or hms!
# LOOKS LIKE HMS!!!
# In which case, looks like the above calculatons are wrong because they assume these are decimal degrees.
# Oh wait, this was already parsed...so it'll look like hms no matter what.

# For debugging purposes, I want to see this organized by aircraft, not by track time
aircraftRecord = {}
for curTime in sorted(trackTimeRecord):
    for curAC in trackTimeRecord[curTime]:
        # print curAC
        callSign = curAC['CallSign']
        if callSign not in aircraftRecord:
            aircraftRecord[callSign] = {'time':[], 'tracks':[], 'plans':[]}

        aircraftRecord[callSign]['time'].append(curTime)
        aircraftRecord[callSign]['tracks'].append(curAC['infoString1'])
        aircraftRecord[callSign]['plans'].append(curAC['infoString2'])

for curAC in aircraftRecord:
    print "{0} -> {1}".format(curAC, aircraftRecord[curAC]['plans'][0])

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

        timeVec         = [obj.text for obj in curTrackPoints.iter('Time')]
        dateVec         = [obj.text for obj in curTrackPoints.iter('Date')]

        datetimeVec     = [datetime.strptime(d + ' ' + t, '%Y-%m-%d %H:%M:%S') for (d,t) in zip(dateVec,timeVec)]
        flightLevelVec  = [obj.text for obj in curTrackPoints.iter('Altitude')]
    except:
        print 'NO TRACKPOINTS -- SKIPPING'
        continue

    # if Callsign not in aircraftRecord:
    #     aircraftRecord[Callsign] = {'time':[], 'tracks':[], 'plans':[]}

    if Callsign in aircraftRecord:
        print "ERROR: {0} has multiple entries".format(Callsign)
        sys.exit()

    aircraftRecord[Callsign] = {'time':timeVec, 'tracks':zip(latVec, lonVec), 'plans':fpVec}






