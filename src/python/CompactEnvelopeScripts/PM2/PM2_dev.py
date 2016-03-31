'''
#
This file creates an envelope for a SpaceX Falcon9 launch
#
'''

'''
System imports that let python know about the file structure of the code
Note: This block should be the same across all main scripts

If you want to use valgrind to debug, compile with -g for debug symbols then call
valgrind --tool=memcheck --suppressions=valgrind-python.py python -E -tt falcon9.py
Be sure to remove -g from compilation when done otherwise code will be slooooow
'''

'''These are the most-likely-to-be-changed parameters'''
freshWind   = True
freshDebris = True
debug       = False

doMain      = True
addStageReentry = False


import os
import sys

# Find the path of the current file, then Point to the root of the package so I can run this script from anywhere
curFilePath = os.path.dirname(os.path.abspath(__file__)) + "/"
rootDir =   os.path.abspath(curFilePath + "../../../../") + "/"
outputDir = rootDir + "outputs/" # Where to store results, gitignored
tempDir =   rootDir + "temp/"   # temp files here, gitignored


'''
Import the modules necessary for the script.
Note: Must Come After sys.path update
Note: This block is probably the same across all main scripts
'''
import CompactEnvelopeBuilder as ceb
from Simulation import TJC
import datetime as dt

from Simulation import LaunchSites
from Simulation import LaunchProviders

from copy import deepcopy


# These parameters get injected into the final footprint name
vehicleName     = LaunchProviders.PM2
launchLocation  = LaunchSites.WSands
vehicleNotes    = 'workInProgress'

# I want to specify the launch location in this way, as opposed to pulling the location from the first state vector,
#   because vehicles that aren't vertical-takeoff may not begin firing until some distance away from the 'launch pad'.
launchLat = LaunchSites.siteDict[launchLocation]['lat']
launchLon = LaunchSites.siteDict[launchLocation]['lon']

# Define statements from inside SkyGrid.h
PROB_IMPACT         = 1001
PROB_CASUALTY       = 1002
PROB_CATASTROPHE    = 1003


'''
# ============= Define the Mission ====================== #
This section is rather convoluted and should be cleaned up.
'''
# from ReadMission import readInput

# Initialize the mission
propagationParamFile = []   # Points to thrust profile for doing propagations
precomputedParamFile = 'PM2_WSands_40kmBasedOnSonda3_COARSENED.txt'                       # Points to file with precomputed profile for nominal trajectory
pathToMissionFiles = './'                       # Kind of a holdover from a previous file structure

# Planet info
omegaE = 7.2921158494529352e-05             # rad/s
planetModel = 0                             # 0 means spherical, 1 means elliptical

curMission = dict(propagationParamFile = propagationParamFile, precomputedParamFile = precomputedParamFile,
                pathToMissionFiles = pathToMissionFiles, omegaE = omegaE, planetModel = planetModel)
curMission = TJC.InitializeMission(curMission)

# If we're propagating trajectories, this needs to be set.  Should do it in the input file, but this is fine for now.
curMission['useLoverD'] = False
curMission['loverd']    = 0.

'''
Defines the file structure for folders / files that are specific to this vehicle.
In theory, all vehicles will have the same file structure and this section would
    be the same across all main files.  Also, should bring pathToMissionFiles
    into this section
'''
# These hold output files, create them in a moment if they don't already exist
curMission['GeneratedFilesFolder']    = curMission['pathToMissionFiles']    + 'GeneratedFiles/'
curMission['debrisPickleFolder']      = curMission['GeneratedFilesFolder']  + 'debrisPickleFolder'
curMission['footprintVectorFolder']   = curMission['GeneratedFilesFolder']  + 'footprintVectorFolder'
curMission['footprintLibrary']        = curMission['GeneratedFilesFolder']  + 'footprintLibrary/'
curMission['facetFolder']             = curMission['GeneratedFilesFolder']  + 'facetFolder/'

# These hold files that need to be read in
curMission['debrisCatPath']           = curMission['pathToMissionFiles'] + 'DebrisCatalog/'
curMission['debrisCatFile']           = 'LynxDebrisCatalog.txt'         # TODO : Custom debris catalog
curMission['atmospherePickle'] = '../AtmoProfiles/WSands.pkl'


# Make sure that the directory for holding the general Generated files exists
folderPath = os.path.abspath(curMission['GeneratedFilesFolder'])
if not os.path.exists(folderPath):
    os.makedirs(folderPath)

# Make sure that debrisPickleFolder exists
folderPath = os.path.abspath(curMission['debrisPickleFolder'])
if not os.path.exists(folderPath):
    os.makedirs(folderPath)

# Make sure that footprintVectorFolder exists
folderPath = os.path.abspath(curMission['footprintVectorFolder'])
if not os.path.exists(folderPath):
    os.makedirs(folderPath)

# Make sure that footprintLibrary exists
folderPath = os.path.abspath(curMission['footprintLibrary'])
if not os.path.exists(folderPath):
    os.makedirs(folderPath)

# Make sure that facetFolder exists
folderPath = os.path.abspath(curMission['facetFolder'])
if not os.path.exists(folderPath):
    os.makedirs(folderPath)

del folderPath      # Just to be safe!




'''
PROPAGATION and PROBABILITY PARAMETERS
Set parameters related to:
    * ASH / density estimation
    * Sky grid granularity
    * NAS reaction time
    * Important time steps (this is kind of confusing)
'''
# Defines the granularity of the gridded sky.  Note, you're restricted to squares in the horizontal plane
# # Parameters for the ASH
NASkm = 18.289

curMission['deltaXY']                   = 0.5    #km
curMission['deltaZ']                    = NASkm/4.   #km
curMission['h1']                        = 3.    # Smoothing parameters for the ASH.  Should be >= deltaXY
curMission['h2']                        = 3.

# Parameters for the safety architecture of the NAS
curMission['reactionTimeMinutes']       = 5     # The number of minutes that the NAS needs to safely handle a sudden debris event.
curMission['thresh']                    = 1e-7  # This is the probability threshold that the cumulative risk must fall below.  Keep in mind
                                                #   there are different definitions of "cumulative" AND there are multiple types of probability.
                                                #   These differences are currently hardcoded and must be changed / recompiled.
curMission['cumulative']                = 'FAA' # The definition for 'cumulative' that we wish to use.
                                                # Options are: FAA, TJC
                                                # Note that if FAA is chosen, the grid will be coarsened to reflect how the FAA calculates
                                                #   hazard areas.  The values of deltaXY and deltaZ will be updated at the appropriate time
                                                #   but they WILL change.

curMission['whichProbability']          = PROB_IMPACT  # Options are IMPACT, CASUALTY, CATASTROPHE

# The different time steps within the mission
curMission['deltaT']                  = 1.      # Seconds, this is the time resolution of a propagated trajectory
curMission['deltaTFail']              = 5.0     # Seconds, this is how often we explode the rocket
curMission['all_points_delta_t']      = 60.0    # Seconds, this will be the time resolution of a compact envelope
                                                #       should be GREATER THAN OR EQUAL to deltaT
curMission['numPiecesPerSample']      = 10      # The number of pieces to consider within each debris group
curMission['useAircraftDensityMap']   = False   # Do we use a uniform or the MIT density map?

curMission['numNodes']                  = 4
curMission['numNodesEnvelopes']         = 4
curMission['NASkm']                     = NASkm

'''
FAILURE PARAMETERS
Import / set parameters related to probabilities of FAILURE for the vehicle
'''
# Generate a realistic profile
from failProfile import failProfile, failProfileSeconds   # This should go in the readInput file
curMission['failProfile'] = failProfile
curMission['failProfileSeconds'] = failProfileSeconds
curMission['pFail'] = 0.02/curMission['all_points_delta_t']     # Probability that vehicle will fail somewhere


'''
Do I still need these?
'''

# Updated, but not sure i need these (this is the first point from the nominal trajecotry file, read this in eventually)
curMission['launchLat'] = launchLat
curMission['launchLon'] = launchLon
curMission['launchAlt'] = 0.   #km

# Not updated, also not sure i need these
curMission['initialUTC'] = 156.84861111111113 # (i think) This number could be anything, as long as it's consistent
# curMission['launchAzimuth'] = 169.    #degrees, this is the heading angle of the SSA runway.  Measured with Google Earth
curMission['launchAzimuth'] = 0.    #Azimuth undefined, basically straight up and down


'''
OUTPUT options
'''
# Export files as GoogleEarth or FACET
curMission['exportGE']    = True
curMission['exportFACET'] = True

# This date gets used for GE.  I believe the FACET files are date agnostic, but i think the time of day might get set here
yyyy    = 2014
mm      = 1
dd      = 1
hour    = 2
min     = 13
sec     = 0
curMission['ExportDate'] = [yyyy, mm, dd, hour, min]
ExportDate = dt.datetime(year=yyyy, month=mm, day=dd, hour=hour, minute=min, second=sec)
curMission['ExportDateDT'] = ExportDate

'''
#
# #
# ======================= Begin Computations ============================ #
# #
#
'''

# Precompute some Atmosphere and Trajectory profiles
freshWind   = False
freshDebris = False
debug       = False

profiles = []
if (freshWind):
    # Should really move all the important mission stuff into this if-statement and wrap it up into the montecarlo dictionary

    numTrajSamples = 1
    numWindSamples = 60

    # I only need to generate wind profiles here, since i'm not going to worry about multiple nominal trajectories yet
    # Could / should probably anticipate doing it though andjust replicate the single trajectory here to conform with the existing infrastrcture

    atmStorage, stateVecStorage, thetagStorage, tfailStorage = TJC.GenerateWindTrajProfiles(curMission, numTrajSamples, numWindSamples)
    profiles = dict(atmStorage = atmStorage, stateVecStorage = stateVecStorage, thetagStorage = thetagStorage, tfailStorage = tfailStorage,
                    numTrajSamples = numTrajSamples, numWindSamples = numWindSamples)

    import pickle
    output = open(curMission['pathToMissionFiles'] + 'localProfiles.pkl', 'wb')
    pickle.dump(profiles,output)
    output.close()

    # tfail = 10.
    # TJC.MonteCarlo_until_tfail(curMission, profiles, tfail)
    # TJC.PlotDebrisFromExplodeTime(curMission, profiles, tfail*1.0)

    # sys.exit()
else:
    import pickle
    profiles = pickle.load(open(curMission['pathToMissionFiles'] + 'localProfiles.pkl','rb'))

if freshDebris:
    t_lo = .0
    t_hi = 230.

    TJC.MonteCarlo_until_tfail(curMission, profiles, t_lo, t_hi)

    ## Find the time until the airspace can become reactive
    #tProactive = TJC.FindStateTimeForProactiveArchitecture(mission1, profiles)
    # sys.exit()


# ## Find the time until the airspace can become reactive
minTime = 30.
maxTime = 90.
tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
print "tProactive = {0}\n".format(tProactive)
# sys.exit()



if not debug:

    # tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles)
    # print "tProactive = {0}\n".format(tProactive)

    footprintStart = 0.
    footprintUntil = 230.
    footprintIntervals = 60.

    footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)
    vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))

    # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
    # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(outfileStr)

    sys.exit()



# if not debug:
#     footprintStart = 0.
#     footprintUntil = 180.
#     footprintIntervals = 60.
#     footprintTotal = []
#
#     for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
#         timelo = footprintStart + ix*footprintIntervals
#         timehi = np.min( (footprintStart + (ix+1)*footprintIntervals - curMission['deltaTFail'], footprintUntil) )
#         EVstrike, Fprint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
#
#         print 'EV =  ' + str(EVstrike)
#
#         Fprint.SmoothedOut(newDeltaT=footprintIntervals)
#         # Fprint.SmoothedOut()
#         Fprint.ExportGoogleEarth('GeneratedFiles/PythonGE_' + str(timelo) + 'To'
#                                               + str(timehi) + 'FootprintSMOOTH.kml', yyyy, mm, dd, hour, min)
#
#         if ix == 0:
#             footprintTotal = Fprint
#         else:
#             footprintTotal.MergeFootprintVectors(Fprint)
#
#
#     # Just to be safe(?), set the params we need in order to translate / rotate
#     footprintTotal.SetAzimuthDeg(curMission['launchAzimuth'])
#     footprintTotal.SetLaunchLatDeg(curMission['launchLat'])
#     footprintTotal.SetLaunchLonDeg(curMission['launchLon'])
#
#
#     vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
#
#     footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
#     outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
#     footprintTotal.StoreFootprintAsVector(outfileStr)
#

