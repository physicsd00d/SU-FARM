'''
#
This file creates an envelope for a Reentry vehicle
#
'''

'''
System imports that let python know about the file structure of the code
Note: This block should be the same across all main scripts

If you want to use valgrind to debug, compile with -g for debug symbols then call
valgrind --tool=memcheck --suppressions=valgrind-python.py python -E -tt falcon9.py
Be sure to remove -g from compilation when done otherwise code will be slooooow
'''
freshMain               = False  # State Vector
freshWind               = True  # Why uncertain wind for this case? B/c uncertainty in direction is manually tweaked.
freshDebris             = True
doMain                  = True

import os
import sys

# Want to import some things that are general to all missions
curFilePath = os.path.dirname(os.path.abspath(__file__)) + "/"      # Reference everything from the location of the file
sys.path.insert(0, os.path.abspath(curFilePath+'../'))              # Back up one so we can import CommonThemes
import CommonThemes as ct

# Find the path of the current file, then Point to the root of the package so I can run this script from anywhere
rootDir =   os.path.abspath(curFilePath + "../../../../") + "/"
outputDir = rootDir + "outputs/" # Where to store results, gitignored
tempDir =   rootDir + ct.tempFolderName   # temp files here, gitignored
debrisPath = rootDir + "src/python/packages/DebrisCatalogs/"


'''
Import the modules necessary for the script.
Note: Must Come After sys.path update
Note: This block is probably the same across all main scripts
'''
# import CompactEnvelopeBuilder as ceb
# import numpy as np
from Simulation import TJC
import datetime as dt

# import matplotlib
# matplotlib.use('Agg')  # Allows plot generation on server without X-windows
# import matplotlib.pyplot as plt

from Simulation import LaunchSites
from Simulation import LaunchProviders

# These parameters get injected into the final footprint name
vehicleName     = LaunchProviders.Reentry
launchLocation  = LaunchSites.OK    # NOTE: even though it says 'launch', in this context it really means 'landing'
vehicleNotes    = ct.vehicleNotes

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
from TrajectoryInfo import initVec_ReentryToOKwithLDpt3 as initVec

# Initialize the mission
propagationParamFile = []                   # Points to thrust profile for doing propagations
precomputedParamFile = 'ReentryToOKwithLDpt3.txt'  # Points to file with precomputed profile for nominal trajectory
pathToMissionFiles = curFilePath            # Kind of a holdover from a previous file structure

# Planet info
omegaE = 7.2921158494529352e-05             # rad/s
planetModel = initVec.planetModel           # 0 means spherical, 1 means elliptical

curMission = dict(propagationParamFile = propagationParamFile, precomputedParamFile = precomputedParamFile,
                pathToMissionFiles = pathToMissionFiles, omegaE = omegaE, planetModel = planetModel)
curMission = TJC.InitializeMission(curMission)
TJC.SetupOutputFolders(curMission, tempDir, outputDir, vehicleName, launchLocation)

# Use the same parameters that were used in the initVec that generated the main trajectories
if initVec.cloption == 0:
    curMission['useLoverD'] = True
else:
    curMission['useLoverD'] = False
curMission['loverd'] = initVec.loverd

# These hold files that need to be read in
curMission['debrisCatPath']           = debrisPath + 'Columbia/'
# curMission['debrisCatFile']           = 'columbiaWithBlast.txt'
curMission['debrisCatFile']           = 'testFileDistributed.txt'
#curMission['debrisCatFile']           = 'debugDistributed.txt'
curMission['atmospherePickle']        = rootDir + "data/AtmoProfiles/WestTexas.pkl"


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

curMission['deltaXY']                   = 2.    #km
curMission['deltaZ']                    = NASkm/1.   #km
curMission['h1']                        = 6.    # Smoothing parameters for the ASH.  Should be >= deltaXY
curMission['h2']                        = 6.

# Parameters for the safety architecture of the NAS
curMission['reactionTimeSeconds']       = 5*60. # The number of seconds that the NAS needs to safely handle a sudden debris event.
                                                #   Negative means to turn off reaction time and keep all points
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
                                                #   Make sure this matches the timestep of the trajectory you have
curMission['deltaTFail']              = 1.     # Seconds, this is how often we explode the rocket
curMission['all_points_delta_t']      = 60.0     # Seconds, this will be the time resolution of a compact envelope
                                                #       should be GREATER THAN OR EQUAL to deltaT
                                                #       For reentry, appears to control the deltaT of the movies made
curMission['numPiecesPerSample']      = 10      # The number of pieces to consider within each debris group
                                                #       IF EMPTY, that means use the actual number for each debris group
curMission['useAircraftDensityMap']   = False   # Do we use a uniform or the MIT density map?
curMission['debrisTimeLimitSec']      = 1*3600  # This is how long to propagate a trajectory for.  If it hasn't landed yet, then give up.
curMission['healthMonitoringLatency'] = 0.      # Seconds

curMission['numNodes']                  = 10  
curMission['numNodesEnvelopes']         = 10
curMission['NASkm']                     = NASkm


'''
FAILURE PARAMETERS
Import / set parameters related to probabilities of FAILURE for the vehicle
'''
# Generate a realistic profile
from failProfile import failProfile, failProfileSeconds   # This should go in the readInput file
curMission['failProfile'] = failProfile
curMission['failProfileSeconds'] = failProfileSeconds
curMission['pFail'] = 0.02    # Probability that vehicle will fail somewhere



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
curMission['launchAzimuth'] = 328.25    #degrees, this is what it looks like on Google Earth
                                        # NOTE: Azimuth here is reversed from launch.
                                        #   In other words, pretend this is actually a launch from the specified site.


'''
OUTPUT options
'''
# Export files as GoogleEarth or FACET
curMission['exportGE']    = False
curMission['exportFACET'] = False

# This date gets used for GE.  I believe the FACET files are date agnostic, but i think the time of day might get set here
yyyy = 2003
mm = 2
dd = 1
hour = 13
min = 59        # min = 49
sec = 30        # sec = 15
ExportDate = dt.datetime(year=yyyy, month=mm, day=dd, hour=hour, minute=min, second=sec)
curMission['ExportDate'] = [yyyy, mm, dd, hour, min]
curMission['ExportDateDT'] = ExportDate

# secondsFromMidnightUTC = hour*3600 + min*60 + sec


'''Unique to Columbia?'''
curMission['isReentry'] = True
# Wind angle measured from East, positive is counterclock
curMission['noWind'] = False
# angleLow    = 0.
# angleHi     = 0.
angleHi         = 360.
angleLow        = 0. #-90
windMagCoeff    = 1.0   # Scale the wind profile


'''
#
# #
# ======================= Begin Computations ============================ #
# #
#
'''

# Can't really incorporate this because to create curMission, the trajectory.txt must already exist
# If you have changed anything in initVec, then you will want to regenerate the nominal trajectory
if (freshMain):
    # import subprocess
    # # subprocess.Popen("Sandbox/GenerateMainPieceTrajectory.py", shell=True)
    # subprocess.call("cd Sandbox && python GenerateMainPieceTrajectory.py && cd -", shell=True)

    from TrajectoryInfo import GenerateMainPieceTrajectory as GMPT
    GMPT.Generate(curMission)
    print "You generated a new trajectory.  Now must restart the script..."
    sys.exit()


curMission['numTrajSamples'] = 1
curMission['numWindSamples'] = 3   # Best results if this is a multiple of the number of nodes you're running on.

profiles = []
if (freshWind):
    # Should really move all the important mission stuff into this if-statement and wrap it up into the montecarlo dictionary

    # I only need to generate wind profiles here, since i'm not going to worry about multiple nominal trajectories yet
    # Could / should probably anticipate doing it though andjust replicate the single trajectory here to conform with the existing infrastrcture

    atmStorage, stateVecStorage, thetagStorage, tfailStorage = \
                            TJC.GenerateWindTrajProfiles(curMission, curMission['numTrajSamples'], curMission['numWindSamples'])
    profiles = dict(atmStorage = atmStorage, stateVecStorage = stateVecStorage, thetagStorage = thetagStorage, tfailStorage = tfailStorage,
                    numTrajSamples = curMission['numTrajSamples'], numWindSamples = curMission['numWindSamples'])

    import pickle
    output = open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl', 'wb')
    pickle.dump(profiles,output)
    output.close()

else:
    import pickle
    profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))


if freshDebris:
    # t_lo = 0.
    # t_hi = 180.
    t_lo = profiles['tfailStorage'][0][0][0] # By setting these to 0 and 1, we'll explode at just the lower time
    t_hi = profiles['tfailStorage'][0][0][-1] # IF YOU CHANGE THIS!!!  Then you'll need to fix the risk calculations later that assume zero
    t_hi = 440.  # IF YOU CHANGE THIS!!!  Then you'll need to fix the risk calculations later that assume zero
#                          # If you set them both to [], then will explode at all times
    TJC.MonteCarloDebris(curMission, profiles, t_lo, t_hi)


# ## Find the time until the airspace can become reactive
#minTime = 0.
#maxTime = 180.
#tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
#print "tProactive = {0}\n".format(tProactive)
# TJC.PlotDebrisFromExplodeTime(curMission, profiles, maxTime, cutoffNAS = False)
# TJC.PlotNominalTrajectories(profiles, curMission, maxTime)
#sys.exit()

footprintIntervals = curMission['all_points_delta_t']
vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))
vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
mainFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '.dat'
totalFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '_stageDown.dat'

if doMain:
    # Note: this is my new and improved method
    curMission['armLength'] = 10000.

    footprintStart = profiles['tfailStorage'][0][0][0] # By setting these to 0 and 1, we'll explode at just the lower time
    # footprintUntil = profiles['tfailStorage'][0][0][-1] # IF YOU CHANGE THIS!!!  Then you'll need to fix the risk calculations later that assume zero
    footprintUntil = 440.

    footprintTotal = TJC.GenerateCompactEnvelopes(curMission, footprintStart, footprintUntil)
    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    footprintTotal.StoreFootprintAsVector(mainFootprintFile)





















# profiles = []
# if (freshWind):
#     # Should really move all the important mission stuff into this if-statement and wrap it up into the montecarlo dictionary
    
#     numTrajSamples = 1      # If you change this to anything other than 1, it might break.  Look at numDebrisPerIXSimulated to start.
#     numWindSamples = 3      # You can increase this, but it will make envelopes smaller.  3 is good enough to replicate paper.

#     # I only need to generate wind profiles here, since i'm not going to worry about multiple nominal trajectories yet
#     # Could / should probably anticipate doing it though andjust replicate the single trajectory here to conform with the existing infrastrcture

#     atmStorage, stateVecStorage, thetagStorage, tfailStorage\
#         = TJC.GenerateWindTrajProfilesDirectional(curMission, numTrajSamples, numWindSamples, angleLow, angleHi, windMagCoeff)

#     profiles = dict(atmStorage = atmStorage, stateVecStorage = stateVecStorage, thetagStorage = thetagStorage,
#                     tfailStorage = tfailStorage, numTrajSamples = numTrajSamples, numWindSamples = numWindSamples)

#     import pickle
#     output = open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl', 'wb')
#     pickle.dump(profiles,output)
#     output.close()

#     # sys.exit()
# else:
#     import pickle
#     profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))


# if freshDebris:
#     # Generate the debris
#     coeffIX = []
#     lowerBreakLimit = profiles['tfailStorage'][0][0][0] # By setting these to 0 and 1, we'll explode at just the lower time
#     upperBreakLimit = profiles['tfailStorage'][0][0][-1] # IF YOU CHANGE THIS!!!  Then you'll need to fix the risk calculations later that assume zero
#                          # If you set them both to [], then will explode at all times
#     TJC.MonteCarlo_Distributed_Reentry_Wrapper_CAIB(curMission, coeffIX, curMission['numPiecesPerSample'],
#                                                     lowerBreakLimit, upperBreakLimit, profiles)

# # ## Find the time until the airspace can become reactive
# #minTime = 0.
# #maxTime = 180.
# #tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
# #print "tProactive = {0}\n".format(tProactive)
# #TJC.PlotNominalTrajectories(profiles, curMission, maxTime)
# #sys.exit()

# footprintIntervals = curMission['all_points_delta_t']
# vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))
# vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
# mainFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '.dat'
# totalFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '_stageDown.dat'

# if doMain:

#     curMission['numWindProfiles'] = len(profiles['atmStorage'])  # Reentry files need this curMission entry
#     # tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles)
#     # print "tProactive = {0}\n".format(tProactive)

#     footprintStart = profiles['tfailStorage'][0][0][0]
#     footprintUntil = profiles['tfailStorage'][0][0][-1] # Landed at this point.
#     print 'footprintStart = {0}'.format(footprintStart)
#     print 'footprintUntil = {0}'.format(footprintUntil)

#     footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)

#     # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
#     # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

#     footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
#     # outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
#     footprintTotal.StoreFootprintAsVector(mainFootprintFile)




# ## Debug
# import pickle
# inFileName = '/Volumes/Storage/Research/SU-FARM/temp/Reentry_OK/debrisPickleFolder/mpc_0_419.pkl'
# input = open(inFileName, 'rb')
# good_mpc = pickle.load(input)
# input.close()

# inFileName = '/Volumes/Storage/Research/SU-FARM/temp/Reentry_OK/debrisPickleFolder/mpc_1_419.pkl'
# input = open(inFileName, 'rb')
# bad_mpc = pickle.load(input)
# input.close()

# curID = 3
# fpArray = bad_mpc['flatPointArray']
# numTimeSteps = bad_mpc['numTimeSteps']
# elemLen = bad_mpc['sizeFlatPointArray']

# fpArray[numTimeSteps.cumsum()[curID-1]*elemLen:].reshape((numTimeSteps[3],elemLen))
# # It landed...must not be properly handling landed points?
# # * Check what pointcloud does with negative-z points
# # * Check what skygrid does
# # * Or could be problem with how I handle the last tstep, not a negative-z problem.

# # Found it.  I have to update the maxtimesteps for each curID to reflect the maximum across all previously passed in pointclouds that had the same curID.

