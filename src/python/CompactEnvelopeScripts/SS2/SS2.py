'''
#
This file creates an envelope for a Lynx MarkII launch and reentry
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
# addStageReentry = False


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
import CompactEnvelopeBuilder as ceb
from Simulation import TJC
import datetime as dt

from Simulation import LaunchSites
from Simulation import LaunchProviders

# from copy import deepcopy


# These parameters get injected into the final footprint name
vehicleName     = LaunchProviders.SS2
launchLocation  = LaunchSites.America
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
# If you do a propagation, then you need to worry about dtval!
propagationParamFile = []                   # Points to thrust profile for doing propagations
precomputedParamFile = 'HTHLwCarrier_Abridged.txt'  # Points to file with precomputed profile for nominal trajectory
pathToMissionFiles = curFilePath            # Kind of a holdover from a previous file structure

# Planet info
omegaE = 7.2921158494529352e-05             # rad/s
planetModel = 0                             # 0 means spherical, 1 means elliptical

curMission = dict(propagationParamFile = propagationParamFile, precomputedParamFile = precomputedParamFile,
                pathToMissionFiles = pathToMissionFiles, omegaE = omegaE, planetModel = planetModel)
curMission = TJC.InitializeMission(curMission)
TJC.SetupOutputFolders(curMission, tempDir, outputDir, vehicleName, launchLocation)


# If we're propagating trajectories, this needs to be set.  Should do it in the input file, but this is fine for now.
curMission['useLoverD'] = False
curMission['loverd']    = 0.


# # These hold files that need to be read in
curMission['debrisCatPath']           = debrisPath + 'Lynx/'
curMission['debrisCatFile']           = 'LynxDebrisCatalog.txt'
curMission['atmospherePickle']  = rootDir + "data/AtmoProfiles/SpaceportAmerica.pkl"



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

curMission['deltaXY']                   = .5    #km
curMission['deltaZ']                    = NASkm/4.   #km
curMission['h1']                        = 3.    # Smoothing parameters for the ASH.  Should be >= deltaXY
curMission['h2']                        = 3.

# Parameters for the safety architecture of the NAS
curMission['reactionTimeSeconds']       = 5*60. # The number of seconds that the NAS needs to safely handle a sudden debris event.
curMission['thresh']                    = 1e-7  # This is the probability threshold that the cumulative risk must fall below.  Keep in mind
                                                #   there are different definitions of "cumulative" AND there are multiple types of probability.
                                                #   These differences are currently hardcoded and must be changed / recompiled.
curMission['cumulative']                = 'FAA' # The definition for 'cumulative' that we wish to use.
                                                # Options are: FAA, TJC
curMission['whichProbability']          = PROB_IMPACT  # Options are IMPACT, CASUALTY, CATASTROPHE

# The different time steps within the mission
curMission['deltaT']                  = 1.      # Seconds, this is the time resolution of a propagated trajectory
                                                # NOTE: This might be REQUIRED to be 1, otherwise holes in PointCloud
                                                # Envelope is half the size if =1 vs =5
                                                # Alternatively, might be required to be deltaTFail because must nest.
curMission['deltaTFail']              = 1.0     # Seconds, this is how often we explode the rocket
# IMPORTANT NOTE: When doing instantaneous health monitoring, if you increase deltaTFail you increase the length of latency
#  with the VHM.  Delta_H = 0 means you always know about all previous timesteps, but if your previous timestep is many
#  seconds away, that could be very noticeable uncertainty.  Further, it loads all the probabilty of failure  of the uncalculated
#  failure times into the failures we did calculate, which makes each explosion about a factor of deltaTFail more risky.
curMission['all_points_delta_t']      = 5.0    # Seconds, this will be the time resolution of a compact envelope
                                                #       should be GREATER THAN OR EQUAL to deltaT
curMission['numPiecesPerSample']      = 10      # The number of pieces to consider within each debris group
curMission['useAircraftDensityMap']   = False   # Do we use a uniform or the MIT density map?
curMission['debrisTimeLimitSec']      = 1*3600  # This is how long to propagate a trajectory for.  If it hasn't landed yet, then give up.
curMission['healthMonitoringLatency'] = 0.      # Seconds

curMission['numNodes']                  = 10 # Will need to install pp to use more nodes
curMission['numNodesEnvelopes']         = 10
curMission['NASkm']                     = NASkm


if curMission['deltaT'] != 1.0:
    print "ERROR: Required deltaT = 1."
    sys.exit()

if (curMission['healthMonitoringLatency'] % curMission['deltaTFail']) != 0.:
    print "ERROR: If you're not exploding every second, then your VHM latency must be a multiple of deltaTFail"
    sys.exit()


'''
FAILURE PARAMETERS
Import / set parameters related to probabilities of FAILURE for the vehicle
'''
# Generate a realistic profile
from failProfile import failProfile, failProfileSeconds   # This should go in the readInput file
curMission['failProfile'] = failProfile
curMission['failProfileSeconds'] = failProfileSeconds
curMission['pFail'] = 0.02     # Probability that vehicle will fail somewhere


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
# curMission['launchAzimuth'] = 180.    #degrees, this is what it looks like on Google Earth
curMission['launchAzimuth'] = 0.    #degrees, Caveat: we're starting the launch at the point where SS2 is dropped
                                    #   so it's actually traveling BACK to the pad (180 degrees) but if I make this
                                    #   0 degrees, then I can rotate the envelope about the launch pad and get the angles
                                    #   to look right.  It's kind of a hack, but only kinda.


#
# curFileName = '/Users/marian/Documents/Research/Prop3Dof/CythonFiles/OtherPythonFiles/FootprintLibrary/workInProgressTJC01/SS2_America_workInProgressHEALTH.dat'
# curFootprint = ceb.PyFootprint(footprintFileName=curFileName)
# curFootprint.SetAzimuthDeg(0.0)
# curFootprint.StoreFootprintAsVector(curFileName)
# sys.exit()




'''
OUTPUT options
'''
# Export files as GoogleEarth or FACET
curMission['exportGE']    = False
curMission['exportFACET'] = False

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

curMission['numTrajSamples'] = 1
curMission['numWindSamples'] = 30   # Best results if this is a multiple of the number of nodes you're running on.

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
    t_lo = .0
    t_hi = 420.

    TJC.MonteCarloDebris(curMission, profiles, t_lo, t_hi)

# # ## Find the time until the airspace can become reactive
# minTime = 60.
# maxTime = 180.
# tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
# print "tProactive = {0}\n".format(tProactive)
# TJC.PlotNominalTrajectories(profiles, curMission, maxTime)
# sys.exit()

footprintIntervals = curMission['all_points_delta_t']
vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))
vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
mainFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '.dat'
totalFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '_stageDown.dat'

if doMain:
    # Note: this is my new and improved method
    curMission['armLength'] = 10000.

    footprintStart = 0.
    footprintUntil = 420.
    footprintTotal = TJC.GenerateCompactEnvelopes(curMission, footprintStart, footprintUntil)
    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    footprintTotal.StoreFootprintAsVector(mainFootprintFile)


