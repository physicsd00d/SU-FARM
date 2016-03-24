# Read in a scenario file and pump out the FACET folders and KML files
import os
# import sys

# Find the path of the current file, then Point to the root of the package so I can run this script from anywhere
curFilePath = os.path.dirname(os.path.abspath(__file__)) + "/"
rootDir =   os.path.abspath(curFilePath + "../../../../") + "/"
outputDir = rootDir + "outputs/" # Where to store results, gitignored
tempDir =   rootDir + "temp/"   # temp files here, gitignored

footprintLibrary    = outputDir + "footprintLibrary/"
KMLFolder           = outputDir + "KMLs/"
FacetFolder         = outputDir + "FacetFiles/"
offsetSeconds = 0.     # Use this if you want to tweak the timing just a little bit.  

# from Simulation import LaunchSites
# from Simulation import LaunchProviders

from Simulation import LaunchSiteLocations as Generator
from CompactEnvelopeBuilder import PyFootprint  #Error expected here for some reason, disregard for now

from datetime import datetime

# scenarioFileFolder = "ScenarioFiles/Sandbox/"
scenarioFileFolder = "ScenarioFiles/SingleDaysUTC/"
# scenarioNameList =   ['2018Low_SingleDay.txt']
# scenarioNameList = ['2018Med_SingleDay.txt']
# scenarioNameList = ['2018High_SingleDay.txt']
# scenarioNameList = ['2025Low_SingleDay.txt']
# scenarioNameList = ['2025Med_SingleDay.txt']
# scenarioNameList = ['2025High_SingleDay.txt']
scenarioNameList = ['2018Low_SingleDay.txt','2018Med_SingleDay.txt','2018High_SingleDay.txt','2025Low_SingleDay.txt','2025Med_SingleDay.txt','2025High_SingleDay.txt']


# scenarioNameList = ['falcon9Real.txt']
# scenarioNameList = os.listdir(scenarioFileFolder)

'''
Load the footprints from the library
'''
footprintDict = dict()
# list all the files in the folder
availableFootprints = os.listdir(footprintLibrary)
for line in availableFootprints:
    # First separate out the kmls from the dats
    [baseName, extension] = line.split('.')
    if (extension.upper() == 'DAT'):
        # This is a dat file (ignore kmls)
        [vehicleName, launchLocation, notes] = baseName.split('_')

        if not footprintDict.has_key(vehicleName):
            # We haven't seen this vehicle yet, so initialize the array
            footprintDict[vehicleName] = [line]
        else:
            # We've seen this one before.  Append.
            footprintDict[vehicleName].append(line)


'''
Run through all the scenario files and create the facet / GE results
'''
for elem in scenarioNameList:
    [scenarioName, ext] = elem.split('.')
    if ext.upper() == 'TXT':

        # scenarioName        = "2018Med_SingleDay"
        # scenarioName        = "2025High_SingleDay"
        scenarioFileName    = scenarioFileFolder + scenarioName + ".txt"

        # Make sure that the directory for holding the facet files exists
        folderPath = os.path.abspath(FacetFolder + scenarioName)
        if not os.path.exists(folderPath):
            os.makedirs(folderPath)

        # Make sure that the directory for holding the KML files exists
        folderPath = os.path.abspath(KMLFolder + scenarioName)
        if not os.path.exists(folderPath):
            os.makedirs(folderPath)

        missionStorage = []
        missionNameStringStorage = ""

        # try:
        inputFile = open(scenarioFileName, 'r')
        for line in inputFile:
            key = line.split()
            if len(key) > 0:  # Make sure line's not empty
                if key[0][0] != "#":  # Ignore the comments
        #             print key

        #             # Construct the datetime object
        #             curDateTime = datetime.strptime(key[0] + "-" + key[1], "%Y-%m-%d-%H:%M")
        #             print curDateTime

                    # Store the current mission data structure
                    curMission = Generator.Mission(launchDateTime = datetime.strptime(key[0] + "-" + key[1], "%Y-%m-%d-%H:%M"),
                                         timeZone     = key[2], vehicle     = key[3],
                                         launchSite   = key[4], azimuth     = key[5],
                                         delayAllowed = key[6])

                    # Store the mission in storage
                    missionStorage.append(curMission)

        # Sort all the launches by time
        missionStorage.sort(key=lambda r: r.launchDateTime);

        # Find the earliest time and use that to define the reference day
        midnight = missionStorage[0].launchDateTime.replace(hour=0, minute=0, second=0)


        # Now that the footprintDict has been created, need to find which file is most appropriate for each curMission
        # curMission = missionStorage[0]

        for curMission in missionStorage:
            '''
            # Prototyping new getMostSimilarFootprint
            '''
            datName = []
            oldLocation = []
            if footprintDict.has_key(curMission.vehicle):
                # Do your magic

                # The available footprints for this vehicle are
                vehiclePrints = footprintDict[curMission.vehicle]

                # If only one file is present, then just use that
                if (len(vehiclePrints) == 1):
                    datName = vehiclePrints[0]
                    oldLocation = datName.split('_')[1]
                    print '\n{0} launching from {1}'.format(curMission.vehicle, oldLocation)

                else:
                    print 'Sort through multiple prints...NET YOT INPLEMINTED!!!'
                    raise RuntimeError


            else:
                print 'ERROR!  You need have at least one .dat file for ever vehicle present'
                print curMission.vehicle
                raise RuntimeError


            # Now have to rotate footprint (if necessary)

            '''
            # Prototyping new getPositionedFootprint
            '''

            curFootprint = PyFootprint(footprintLibrary + datName, asVector=True)
            targetAz    = float(curMission.azimuth)
            currentAz   = curFootprint.GetAzimuthDeg()
            # print '     Azimuth current = {0}, target = {1}'.format(currentAz, targetAz)

            if targetAz > 360.:
                print '     Azimuth    : DEFAULT {0}'.format(currentAz)
                print 'CANNOT USE THIS OPTION UNTIL THE DAT FILES HAVE THE PROPER AZIMUTH RECORDED IN THEM'
            else:
                curFootprint.ChangeAzimuthToDeg(float(targetAz))
                print '     Azimuth    : Rotating from {0} to {1}'.format(currentAz, targetAz)
                # raise RuntimeError

            # Set to lat/lon of desired launch site
            targetLocation = curMission.launchSite
            if oldLocation != targetLocation:
                [lat, lon] = Generator.getLatLon(curMission)
                oldLat = curFootprint.GetLaunchLatDeg()
                oldLon = curFootprint.GetLaunchLonDeg()

                # print 'Translating from ' + str([oldLat, oldLon]) + ' to ' + str([lat,lon])
                print '     LaunchSite : Translating from {0} to {1}'.format(oldLocation, targetLocation)
                curFootprint.ChangeLaunchSiteToDeg(lat,lon)
            else:
                print '     LaunchSite : DEFAULT {0}'.format(oldLocation)


            now = curMission.launchDateTime

            # Generate the standard name
            missionDateTime = curMission.launchDateTime.strftime("%Y-%m-%d-%H-%M")
            missionName = curMission.vehicle + "_" + curMission.launchSite + '_' + missionDateTime

            # Create the GE file so we can visually inspect
            GEFileName = str(KMLFolder + scenarioName + '/' + missionName + '.kml')

            # Make the FACET files for this mission
            facetFolderName = str(FacetFolder + scenarioName + "/" + missionName + "/")

            obsolete = -1.
            # minutessFromMidnight = now.hour*60. + now.minute
            # curFootprint.MakeFacetFiles(facetFolderName,minutessFromMidnight,offsetSeconds,obsolete)

            # secondsFromMidnight = now.hour*3600. + now.minute*60.
            secondsFromMidnight = int((now - midnight).total_seconds())
            print "secondsFromMidnight = {0}".format(secondsFromMidnight)
            curFootprint.MakeFacetFiles(facetFolderName,secondsFromMidnight,offsetSeconds,obsolete)

            curFootprint.ExportGoogleEarth(GEFileName, now.year , now.month, now.day, now.hour, now.minute)

print "Done"






















