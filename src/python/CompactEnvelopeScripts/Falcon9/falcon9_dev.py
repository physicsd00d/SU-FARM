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
freshWind   = False
freshDebris = False
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
vehicleName     = LaunchProviders.Falcon9
launchLocation  = LaunchSites.Cape
vehicleNotes    = ''

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
propagationParamFile = 'nominalParam_new.txt'   # Points to thrust profile for doing propagations
precomputedParamFile = []                       # Points to file with precomputed profile for nominal trajectory
pathToMissionFiles = curFilePath                # Kind of a holdover from a previous file structure

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


# These hold files that need to be read in
curMission['debrisCatPath']     = curMission['pathToMissionFiles'] + 'DebrisCatalog/'
curMission['debrisCatFile']     = 'Halcon9_1stNEW.txt'
curMission['atmospherePickle']  = rootDir + "data/AtmoProfiles/Cape.pkl"



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
curMission['reactionTimeSeconds']       = 5*60.     # The number of seconds that the NAS needs to safely handle a sudden debris event.
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
curMission['deltaTFail']              = 10.0     # Seconds, this is how often we explode the rocket
# IMPORTANT NOTE: When doing instantaneous health monitoring, if you increase deltaTFail you increase the length of latency
#  with the VHM.  Delta_H = 0 means you always know about all previous timesteps, but if your previous timestep is many
#  seconds away, that could be very noticeable uncertainty.  Further, it loads all the probabilty of failure  of the uncalculated
#  failure times into the failures we did calculate, which makes each explosion about a factor of deltaTFail more risky.
curMission['all_points_delta_t']      = 60.0    # Seconds, this will be the time resolution of a compact envelope
                                                #       should be GREATER THAN OR EQUAL to deltaT
curMission['numPiecesPerSample']      = 10      # The number of pieces to consider within each debris group
curMission['useAircraftDensityMap']   = False   # Do we use a uniform or the MIT density map?
curMission['debrisTimeLimitSec']      = 1*3600  # This is how long to propagate a trajectory for.  If it hasn't landed yet, then give up.
curMission['healthMonitoringLatency'] = 0.      # Seconds

curMission['numNodes']                  = 8 # Will need to install pp to use more nodes
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
curMission['launchAzimuth'] = 47.8      #degrees, this is what it looks like on Google Earth

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

if debug:
    # Change a few values
    curMission['debrisCatFile']           = 'Debug.txt'
    curMission['reactionTimeSeconds']       = 5*60.     # The number of minutes that the NAS needs to safely handle a sudden debris event.
    curMission['numPiecesPerSample']      = 2      # The number of pieces to consider within each debris group


profiles = []
if (freshWind):
    # Should really move all the important mission stuff into this if-statement and wrap it up into the montecarlo dictionary

    numTrajSamples = 1
    numWindSamples = 2

    # I only need to generate wind profiles here, since i'm not going to worry about multiple nominal trajectories yet
    # Could / should probably anticipate doing it though andjust replicate the single trajectory here to conform with the existing infrastrcture

    atmStorage, stateVecStorage, thetagStorage, tfailStorage = TJC.GenerateWindTrajProfiles(curMission, numTrajSamples, numWindSamples)
    profiles = dict(atmStorage = atmStorage, stateVecStorage = stateVecStorage, thetagStorage = thetagStorage, tfailStorage = tfailStorage,
                    numTrajSamples = numTrajSamples, numWindSamples = numWindSamples)

    import pickle
    output = open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl', 'wb')
    pickle.dump(profiles,output)
    output.close()

    # tfail = 10.
    # TJC.MonteCarlo_until_tfail(curMission, profiles, tfail)
    # TJC.PlotDebrisFromExplodeTime(curMission, profiles, tfail*1.0)

    # sys.exit()
else:
    import pickle
    profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))


if freshDebris:
    t_lo = 0.
    t_hi = 180.

    TJC.MonteCarlo_until_tfail(curMission, profiles, t_lo, t_hi)

# # ## Find the time until the airspace can become reactive
# minTime = 120.
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

# import numpy as np
# t_lo = .0
# t_hi = 170.
# deltaTFail = curMission['deltaTFail']
# timeVec = np.arange(t_hi*1.0,t_lo-deltaTFail,-deltaTFail)        #curTime is in seconds
# for curTime in timeVec:
#     print curTime
#     # TJC.PlotDebrisFromExplodeTime(curMission, profiles, curTime, cutoffNAS = True)
#     TJC.PlotSubEnvelopes(curMission, curTime)
# sys.exit()



# debugSingleTime = False
# if debugSingleTime:
#      ### ===== DEBUG =========
#     tfailSec = 70.

#     from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud#, PyFootprint
#     import pickle
#     import numpy as np

#     deltaXY                 = curMission['deltaXY']
#     deltaZ                  = curMission['deltaZ']
#     h1                      = curMission['h1']
#     h2                      = curMission['h2']
#     debrisPickleFolder      = curMission['debrisPickleFolder']
#     footprintVectorFolder   = curMission['footprintVectorFolder']
#     thresh                  = curMission['thresh']
#     cumulative              = curMission['cumulative']
#     whichProbability        = curMission['whichProbability']

#     inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
#     input = open(inFileName, 'rb')
#     cur_mpc = pickle.load(input)
#     input.close()

#     arefMeanList = cur_mpc['arefMeanList']
#     numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']

#     # Package them up into a PointCLoud
#     # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
#     curPointCloud = PyPointCloud(cur_mpc, tfailSec, curMission)

#     # Place the cloud into a Grid
#     curSkyGrid    = PySkyGrid(curMission=curMission, pointCloud=curPointCloud)

#     # # Now if I ASH without ASHing, that should just give me the unspread probabilities
#     print 'ASHING'
#     h1                        = curMission['deltaXY']     # Smoothing parameters for the ASH.  Should be >= deltaXY
#     h2                        = curMission['deltaXY'] 
#     # h1                        = curMission['h1']     # Smoothing parameters for the ASH.  Should be >= deltaXY
#     # h2                        = curMission['h2'] 
#     curSkyGrid.generateASH(h1, h2)

#     def checkNorm(ash):
#         curNorm = 0.
#         for curZ in ash:
#             for curX in ash[curZ]:
#                 for curY in ash[curZ][curX]:
#                     curNorm += ash[curZ][curX][curY]
#         return curNorm

#     # Okay, now I can look through the histograms any way I want
#     # curID = 10  # highest beta
#     curID = 2   # most pieces.  This must SURELY generate a hazard area.  Very light, mostly hangs in air.
#     hist = dict()
#     ash = dict()
#     whichProb = 0   # Impact
#     for tx in range(300):
#         hist[tx] = curSkyGrid.SendHistogramToPython(curID,tx)
#         ash[tx] = curSkyGrid.SendProbabilitiesToPython(curID,tx, 0)

#         # if len(hist[tx]) > 0:
#         # print "{0}: {1} --> {2}".format(tx, hist[tx], ash[tx])
#         print "{0}: {1} --> {2}".format(tx, hist[tx], 1-checkNorm(ash[tx]))





# tfailSec = 100.
# inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
# input = open(inFileName, 'rb')
# cur_mpc = pickle.load(input)
# input.close()

# TJC.PlotDebrisFromExplodeTime(curMission, profiles, tfail=100., cutoffNAS = False)


if doMain:
    # Note: this is my new and improved method
    curMission['armLength'] = 10000.

    footprintStart = 0.
    footprintUntil = min(180., curMission['failProfileSeconds'][-1])
    footprintTotal = TJC.GenerateCompactEnvelopes(curMission, footprintStart, footprintUntil)
    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    footprintTotal.StoreFootprintAsVector(mainFootprintFile)



if addStageReentry:
    # Note: this is being done with the old method
    '''# Prototype handling the first stage reentry'''

    # Specify the time of the staging
    tStage = 175

    curMission['debrisPickleFolder']      = curMission['GeneratedFilesFolder']  + 'debrisPickleFolder'

    firstStageMission = deepcopy(curMission)
    firstStageMission['debrisPickleFolder'] = curMission['GeneratedFilesFolder']  + 'firstStagePickleFolder'
    firstStageMission['debrisCatFile'] = 'firstStage.txt'
    firstStageMission['reactionTimeSeconds'] = -1
    firstStageMission['numPiecesPerSample'] = 1
    firstStageMission['thresh'] = 0.
    coeffIX = []

    mpc = TJC.MonteCarlo_at_tfail(firstStageMission, coeffIX, tStage, firstStageMission['numPiecesPerSample'], profiles)
    debrisPickleFolder = firstStageMission['debrisPickleFolder']

    # print 'COMMENTED OUT WRITING TO FILE FOR DEBUGGING PURPOSES'
    # Make sure that the output directory exists
    folderPath = os.path.abspath(debrisPickleFolder)
    if not os.path.exists(folderPath):
        os.makedirs(folderPath)

    output = open(folderPath + '/mpc_' + str(tStage) + '.pkl', 'wb')
    pickle.dump(mpc,output,2)
    output.close()

    curPFail = 1.
    EV_strike, outfileStr = TJC.genFootprint(firstStageMission, tStage, curPFail)
    firstStageFootprint = ceb.PyFootprint(footprintFileName=outfileStr)

    firstStageFootprint.SmoothedOut()   #I believe this will simply smooth the footprints and not alter the timesteps

    # Just to be safe(?), set the params we need in order to translate / rotate
    firstStageFootprint.SetAzimuthDeg(firstStageMission['launchAzimuth'])
    firstStageFootprint.SetLaunchLatDeg(firstStageMission['launchLat'])
    firstStageFootprint.SetLaunchLonDeg(firstStageMission['launchLon'])

    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, "firstStage")
    outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    firstStageFootprint.StoreFootprintAsVector(outfileStr)

    # ''' Now make a GE animation of the first stage spread'''
    # myTraj = ceb.PyTrajectory()
    # secondsFromLaunch = tStage
    # myTraj.loadDebrisTrajectory(mpc, secondsFromLaunch, firstStageMission)
    #
    # outFileName = curMission['GeneratedFilesFolder'] + "firstStage.kml"
    # myTraj.ExportGoogleEarth(outFileName, ExportDate)


    '''Now Merge with main footprint'''
    totalFootprint = ceb.PyFootprint(footprintFileName=mainFootprintFile)
    totalFootprint.MergeFootprintVectors(firstStageFootprint)
    totalFootprint.StoreFootprintAsVector(totalFootprintFile)
    totalFootprint.ExportGoogleEarth(firstStageMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)

# t = 0.0, curPFail = 0.00384158054763, curEV = 8.74345357448e-11
# t = 10.0, curPFail = 0.00322909977888, curEV = 5.22387577793e-08
# t = 20.0, curPFail = 0.00167566498909, curEV = 9.51254676391e-08
# t = 30.0, curPFail = 0.000763743561604, curEV = 2.40607358433e-10
# t = 40.0, curPFail = 0.000315583572514, curEV = 9.26786832273e-08
# t = 50.0, curPFail = 0.000122854360758, curEV = 7.65799722236e-08
# t = 60.0, curPFail = 7.79467850187e-05, curEV = 9.46252965686e-08
# t = 70.0, curPFail = 0.000203184252099, curEV = 8.03477871841e-08
# t = 80.0, curPFail = 0.000637036973329, curEV = 9.92331713844e-08
# t = 90.0, curPFail = 0.00144568846946, curEV = 9.98893278579e-08
# t = 100.0, curPFail = 0.002300990413, curEV = 9.95909693532e-08
# t = 110.0, curPFail = 0.00253199788052, curEV = 9.98862611894e-08
# t = 120.0, curPFail = 0.00184510768573, curEV = 9.97920312481e-08
# t = 130.0, curPFail = 0.000813289770532, curEV = 9.98483287621e-08
# t = 140.0, curPFail = 0.000181736694585, curEV = 9.9951091907e-08
# t = 150.0, curPFail = 1.43309723121e-05, curEV = 3.28934668997e-08
# t = 160.0, curPFail = 1.63286854877e-07, curEV = 9.18756160035e-13
# t = 170.0, curPFail = 6.08574524108e-12, curEV = 0.0

