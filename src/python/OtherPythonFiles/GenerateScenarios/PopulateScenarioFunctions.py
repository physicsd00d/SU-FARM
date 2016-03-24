from LaunchSiteLocations import Mission
import random
import math
from datetime import datetime
import calendar
import numpy as np

def distributeLaunches(launchProvider, yyyy):
    # Calculate the reference date
    refOrdinal = datetime.today().replace(day=1, month=1, year=yyyy).toordinal()
    
    # For now assume any day is as likely as any other
    numDaysThisYear = calendar.isleap(yyyy)*366 + (not calendar.isleap(yyyy))*365
        
    # Find the vehicle names
    vehicleNames = launchProvider.keys()
    MissionStorage = []

    for curVehicle in vehicleNames:
        print curVehicle
        curMissionStats = launchProvider[curVehicle]
        
        vehicleName = curMissionStats.vehicleName
        numLaunches = curMissionStats.numLaunches
        launchSites = curMissionStats.launchSites
        launchSiteFrequencies = curMissionStats.launchSiteFrequencies
        
        if np.abs(np.sum(launchSiteFrequencies)-1) > 1e-15:
            print 'ERROR!!! YOUR LAUNCH FREQUENCIES DO NOT SUM TO 1 ~~~~~~~~~~~~~~~~~'
            print str(np.sum(launchSiteFrequencies)-1)
            print str(launchSiteFrequencies)
     
        # Which constraints are active?
        daySpacingConstraintIsActive = curMissionStats.constraintDict['daySpacingConstraint']['isActive']
        daySpacingConstraintValue = curMissionStats.constraintDict['daySpacingConstraint']['value']
        

     
     
        numFailed = 0;
        for ix in range(0,numLaunches):
            # Assign the launch site
            randNum = random.random()
            whichIndex = randNum < np.cumsum(launchSiteFrequencies)
#             print str(randNum) + "  " + str(np.cumsum(launchSiteFrequencies)) + "    " + str(whichIndex)
            randomSite = launchSites[whichIndex][0]
             
            # Assign the launch azimuth
            randomAzimuth = random.random()*180 - 90
             
            # Assign the allowable launch delay
            randomDelayMin = random.randint(5,20)   #Random int between A and B
             
            # Assemble the full random date
            attempts = 0
            minDistance = -1;   # set to -1 to allow you inside the while loop
            #         constraint = 7  #days
            allConstraintsMet = False
         
            maxAttempts = 50
            while ((attempts < maxAttempts) and not allConstraintsMet):
                # Keep track of how many times you tried
                attempts = attempts+1
                 
                # Set to true, then AND with all the contraints met.  If any not met, will trip to False and stay false.
                allConstraintsMet = True
                 
                # Randomly place a launch
                randomDayNumber = random.randint(0,numDaysThisYear-1)   #Days since beginning of year
                launchMinFromMidnight = random.randint(0,60*24-1)       #Minutes since midnight of random day
                hh = int(math.floor(launchMinFromMidnight/60))              #turn that into hours
                mm = launchMinFromMidnight%60                               #and remaining minutes
                randomDate = datetime.fromordinal(refOrdinal + randomDayNumber).replace(hour=hh, minute=mm)
                 
                # Check constraint that days aren't too close to each other
                minDistance = 1e8     #Constraint doesn't apply for ix == 0
                 
                if (ix > 0):
                    for testLaunch in MissionStorage:
                        testDate = testLaunch.launchDateTime
                        deltaDays = abs((randomDate - testDate).days)
                        if deltaDays < minDistance:
                            minDistance = deltaDays
                             
                     
                # Check the constraints
                if daySpacingConstraintIsActive:
                    daySpacingConstraintMet = (minDistance >= daySpacingConstraintValue)
                    allConstraintsMet = allConstraintsMet and daySpacingConstraintMet
                 
             
                 
        #         print "minDistance = " + str(minDistance)
        
            # Save the date if you were able to find one that worked
            if (attempts == maxAttempts):
                print str(minDistance) + ":   Failed to find a date that satisfied the constrants\n"
                numFailed = numFailed +1
            else:
                curMission = Mission(launchDateTime = randomDate,
                                     timeZone   = 'GMT',                    vehicle = vehicleName,
                                     launchSite   = randomSite,             azimuth = round(randomAzimuth,2),
                                     delayAllowed = randomDelayMin)
                MissionStorage.append(curMission)
         
         
        # Make sure nothing failed
        if numFailed != 0:
            print "numFailed = " + str(numFailed)
     
    # Sort the launches by date and time
    MissionStorage.sort(key=lambda r: r.launchDateTime)
     
    return MissionStorage




# Print the Specific Mission
def printSingleMission(curMission):
    for cur in curMission:
        outString = (str(cur.launchDateTime.date()) + "\t"  + cur.launchDateTime.strftime("%H:%M") + "\t" + 
                    cur.timeZone + "\t"     + cur.vehicle + "\t\t" + 
                    cur.launchSite + "\t\t"   + str(cur.azimuth) + "\t\t" + 
                    str(cur.delayAllowed))
        print outString
        
        
        
        
        
        
        
def printMissionRecord(missionRecord):
    # First consolidate everything into one array
    allLaunches = []
    for curProvider in missionRecord:
        for curMission in missionRecord[curProvider]:
            allLaunches.append(curMission)
            
    print allLaunches

    # Sort the launches by date and time
    allLaunches.sort(key=lambda r: r.launchDateTime)
    
    print "Number of launches = " + str(len(allLaunches))

    
    totalOut = ""    
    
    for cur in allLaunches:
        outString = (str(cur.launchDateTime.date()) + "\t"  + cur.launchDateTime.strftime("%H:%M") + "\t\t" + 
                    cur.timeZone + "\t\t\t"     + cur.vehicle + "\t\t\t" + 
                    cur.launchSite + "\t\t"   + '%.2f' % cur.azimuth + "\t\t\t" + 
                    str(cur.delayAllowed))
        totalOut += outString + "\n"
#         print outString
#     print totalOut
    return totalOut