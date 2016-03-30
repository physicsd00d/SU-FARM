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

doMain      = False
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
vehicleNotes    = 'space2015'

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
curMission['healthMonitoringLatency'] = 10.      # Seconds

curMission['numNodes']                  = 4 # Will need to install pp to use more nodes
curMission['numNodesEnvelopes']         = 2
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


# ### Prototyping the use of non-zero health monitoring.

# # Reaction time is applied within PointCloud, so all envelopes are automatically for t <= J

# # This makes a footprint for a failure at ix over all time up to the reaction time, f + delta_R
# # genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) 
# # So that's almost the argument of the second sum, but the pointcloud will need to keep up to f + delta_R + delta_H
# # TODO: Change PointCloud to incorporate f + delta_R + delta_H

# delta_H = curMission['healthMonitoringLatency']/curMission['deltaT']    # TODO round to integer
# delta_R = curMission['reactionTimeSeconds']/curMission['deltaT']  # TODO round to integer
# # Flow should go like this.  Knowing delta_R and delta_H ahead of time
# # time = 0
# # SkyGrid for P_I(x, f=0 | t <= f=0 + delta_R + delta_H) to be used immediately
# # SkyGrid for P_I(x, f=0 | t <= f=0 + delta_H) to be used once the health update comes in
# debrisPickleFolder      = curMission['debrisPickleFolder']
# deltaXY                 = curMission['deltaXY']
# deltaZ                  = curMission['deltaZ']
# tfailSec = 60.

# inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
# input = open(inFileName, 'rb')
# cur_mpc = pickle.load(input)
# input.close()

# arefMeanList = cur_mpc['arefMeanList']
# numberOfPiecesMeanList  = cur_mpc['numberOfPiecesMeanList']

# from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud, PyGrid3D#, PyFootprint

# # Package them up into a PointCLoud
# # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
# curPointCloud           = PyPointCloud(cur_mpc, tfailSec, curMission)

# # Place the cloud into a Grid
# curSkyGrid              = PySkyGrid(curMission, curPointCloud)
# h1                        = curMission['deltaXY']     # Smoothing parameters for the ASH.  Should be >= deltaXY
# h2                        = curMission['deltaXY'] 
# # h1                        = curMission['h1']     # Smoothing parameters for the ASH.  Should be >= deltaXY
# # h2                        = curMission['h2'] 
# curSkyGrid.generateASH(h1, h2)

# def checkNorm(ash):
#     curNorm = 0.
#     for curZ in ash:
#         for curX in ash[curZ]:
#             for curY in ash[curZ][curX]:
#                 curNorm += ash[curZ][curX][curY]
#     return curNorm

# # Okay, now I can look through the histograms any way I want
# # curID = 10  # highest beta
# curID = 2   # most pieces.  This must SURELY generate a hazard area.  Very light, mostly hangs in air.
# hist = dict()
# ash = dict()
# whichProb = 0   # Impact
# for tx in range(300):
#     hist[tx] = curSkyGrid.SendHistogramToPython(curID,tx)
#     ash[tx] = curSkyGrid.SendProbabilitiesToPython(curID,tx, 0)

#     # if len(hist[tx]) > 0:
#     # print "{0}: {1} --> {2}".format(tx, hist[tx], ash[tx])
#     print "{0}: {1} --> {2}".format(tx, hist[tx], 1-checkNorm(ash[tx]))

# print 'generateHazardProbabilities'
# curSkyGrid.generateHazardProbabilities(numberOfPiecesMeanList)

# whichProb = curSkyGrid.getProbImpactCode()
# tempGrid3D = curSkyGrid.GenerateSpatialProbability(whichProb, tfailSec+150, tfailSec)
# temp = tempGrid3D.getGrid()

# for z in temp:
#     for x in temp[z]:
#         for y in temp[z][x]:
#             print "[{0}][{1}][{2}] = {3:e}".format(z,x,y, temp[z][x][y]) 


# newGrid = PyGrid3D(tempGrid3D) # Not okay

# addGrid = tempGrid3D + newGrid








### Prototyping the use of non-zero health monitoring.

# Reaction time is applied within PointCloud, so all envelopes are automatically for t <= J

# This makes a footprint for a failure at ix over all time up to the reaction time, f + delta_R
# genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) 
# So that's almost the argument of the second sum, but the pointcloud will need to keep up to f + delta_R + delta_H
# TODO: Change PointCloud to incorporate f + delta_R + delta_H

delta_H = curMission['healthMonitoringLatency']/curMission['deltaT']    # TODO round to integer
delta_R = curMission['reactionTimeSeconds']/curMission['deltaT']  # TODO round to integer

debrisPickleFolder      = curMission['debrisPickleFolder']
# deltaXY                 = curMission['deltaXY']
# deltaZ                  = curMission['deltaZ']
from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud, PyGrid3D#, PyFootprint




# Flow should go like this.  Knowing delta_R and delta_H ahead of time
# time = 0
# SkyGrid for P_I(x, f=0 | t <= f=0 + delta_R + delta_H) to be used immediately
# SkyGrid for P_I(x, f=0 | t <= f=0 + delta_H) to be used once the health update comes in
tfailSec = 0.

inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
input = open(inFileName, 'rb')
cur_mpc = pickle.load(input)
input.close()

# Package them up into a PointCLoud
curPointCloud           = PyPointCloud(cur_mpc, tfailSec, curMission)

# Place the cloud into a Grid
curSkyGrid              = PySkyGrid(curMission, curPointCloud)

# ASH them
h1                        = curMission['deltaXY']     # Smoothing parameters for the ASH.  Should be >= deltaXY
h2                        = curMission['deltaXY'] 
# h1                        = curMission['h1']     # Smoothing parameters for the ASH.  Should be >= deltaXY
# h2                        = curMission['h2'] 
curSkyGrid.generateASH(h1, h2)

# Calculate all of the hazard probabilities
curSkyGrid.generateHazardProbabilities(cur_mpc['numberOfPiecesMeanList'])

# Finally get the probabilities for the times we want
whichProb = curSkyGrid.getProbImpactCode()
P_RH = curSkyGrid.GenerateSpatialProbability(whichProb, tfailSec + delta_R + delta_H, tfailSec)
P_H = curSkyGrid.GenerateSpatialProbability(whichProb, tfailSec + delta_H, tfailSec)


#
# time = 1
# SkyGrid for P_I(x, f=1 | t <= f=1 + delta_R + delta_H) to be used immediately
# SkyGrid for P_I(x, f=1 | t <= f=1 + delta_H) to be used once the health update comes in
#








# footprintTotal = []

# for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
#     timelo = footprintStart + ix*footprintIntervals
#     timehi = np.min( (footprintStart + (ix+1)*footprintIntervals, footprintUntil) )

#     print 'TIMES: From {0} to {1}'.format(timelo, timehi)
#     EVstrike, curFootPrint = makeFootprintFromTimes_InstantaneousOnly(curMission, timelo, timehi)
#     print 'EV =  ' + str(EVstrike)

#     # Now take that footprint and...
#     # Smooth it out to a single timestep

#     # This is what the footprint will look like if there was an "error" in generating it
#     if curFootPrint == []:
#         print "\n\n===========\nEMPTY FOOTPRINT!!!\n==============\n\n"
#         continue

#     numRange = curFootPrint.getNumRange()
#     curFootPrint.SmoothedOut(numRange * curMission['all_points_delta_t'])  # This will make footprintDelaT = numRange, and then change numRange to = 1

#     numRange = curFootPrint.getNumRange()
#     FPDeltaT = curFootPrint.getDeltaT()

#     # Resize the deltaT to be only the length of the interval
#     #   So if we're making an envelope at each second, then the footprint should be chopped at 1 second
#     #   If we're combining times, like every 5 seconds or every minute, then it should be 5s or 60s
#     curFootPrint.ChopTimeAt(footprintIntervals)

#     # Translate the footprint forward to tfailSec
#     if timelo > 0:
#         curFootPrint.SlideFootprintBySeconds(timelo)

#     # Merge it with the others
#     if ix == 0:
#         footprintTotal = curFootPrint
#     else:
#         print '\n\nMERGE'
#         footprintTotal.MergeFootprintVectors(curFootPrint)


#     # # Print to GE
#     # debugFolder = 'GeneratedFiles/Sandbox/'
#     # vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, timelo, timehi)
#     # curFootPrint.ExportGoogleEarth(debugFolder + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)

#     # # Fprint.SmoothedOut(footprintIntervals)
#     #
#     # # Fprint.SmoothedOut()
#     # curFootPrint.ExportGoogleEarth('GeneratedFiles/PythonGE_' + str(timelo) + 'To'
#     #                                       + str(timehi) + 'FootprintSMOOTH.kml', yyyy, mm, dd, hour, min)

# footprintTotal.SmoothedOut(0)   #I believe this will simply smooth the footprints and not alter the timesteps

# # Just to be safe(?), set the params we need in order to translate / rotate
# footprintTotal.SetAzimuthDeg(curMission['launchAzimuth'])
# footprintTotal.SetLaunchLatDeg(curMission['launchLat'])
# footprintTotal.SetLaunchLonDeg(curMission['launchLon'])

# return footprintTotal





### Match this
# generateHazardProbabilities
# [124][1][-14116][6352] -> 0.000000E+00 * 9.999613E-01
# [125][1][-14116][6352] -> 3.873102E-05 * 9.999614E-01
# [177][1][-14116][6352] -> 7.734734E-05 * 9.977185E-01
# [178][1][-14116][6352] -> 2.358637E-03 * 9.954200E-01
# [179][1][-14116][6352] -> 6.927807E-03 * 9.954235E-01
# [180][1][-14116][6352] -> 1.147262E-02 * 9.954267E-01
# [181][1][-14116][6352] -> 1.599340E-02 * 9.931317E-01
# [182][1][-14116][6352] -> 2.275182E-02 * 9.954342E-01
# [183][1][-14116][6352] -> 2.721378E-02 * 9.463568E-01
# [184][1][-14116][6352] -> 7.939714E-02 * 9.308014E-01
# [185][1][-14116][6352] -> 1.431016E-01 * 9.329589E-01
# [186][1][-14116][6352] -> 2.005490E-01 * 8.832070E-01
# [187][1][-14116][6352] -> 2.939193E-01 * 8.593692E-01
# [188][1][-14116][6352] -> 3.932160E-01 * 8.461824E-01
# [189][1][-14116][6352] -> 4.865500E-01 * 7.613635E-01
# [190][1][-14116][6352] -> 6.090779E-01 * 6.300726E-01
# [191][1][-14116][6352] -> 7.536907E-01 * 5.542229E-01
# [192][1][-14116][6352] -> 8.634898E-01 * 5.403972E-01
# [193][1][-14116][6352] -> 9.262302E-01 * 4.919908E-01
# [194][1][-14116][6352] -> 9.637060E-01 * 4.470612E-01
# [195][1][-14116][6352] -> 9.837743E-01 * 3.884570E-01
# [196][1][-14116][6352] -> 9.936970E-01 * 3.572080E-01
# [197][1][-14116][6352] -> 9.977485E-01 * 3.271343E-01
# [198][1][-14116][6352] -> 9.992635E-01 * 2.962286E-01
# [199][1][-14116][6352] -> 9.997818E-01 * 2.765135E-01
# [200][1][-14116][6352] -> 9.999397E-01 * 2.281497E-01
# [201][1][-14116][6352] -> 9.999862E-01 * 1.329828E-01
# [202][1][-14116][6352] -> 9.999982E-01 * 8.110280E-02
# [203][1][-14116][6352] -> 9.999999E-01 * 6.712495E-02
# [204][1][-14116][6352] -> 1.000000E+00 * 6.856211E-02
# [205][1][-14116][6352] -> 1.000000E+00 * 7.275751E-02
# [206][1][-14116][6352] -> 1.000000E+00 * 7.370047E-02
# [207][1][-14116][6352] -> 1.000000E+00 * 7.514612E-02
# [208][1][-14116][6352] -> 1.000000E+00 * 7.866974E-02
# [209][1][-14116][6352] -> 1.000000E+00 * 7.968113E-02
# [210][1][-14116][6352] -> 1.000000E+00 * 8.333504E-02
# [0][-14119][6352] = 7.222048e-02
# [0][-14119][6353] = 3.279850e-02
# [0][-14119][6354] = 5.102704e-04
# [0][-14118][6352] = 9.982917e-01
# [0][-14118][6353] = 5.089480e-01
# [0][-14118][6354] = 1.050376e-02
# [0][-14118][6350] = 1.263690e-01
# [0][-14118][6351] = 8.839138e-01
# [0][-14117][6352] = 9.999909e-01
# [0][-14117][6353] = 8.021458e-01
# [0][-14117][6354] = 5.295458e-02
# [0][-14117][6351] = 6.394572e-01
# [0][-14116][6352] = 9.999549e-01
# [0][-14116][6353] = 7.140744e-01
# [0][-14116][6354] = 2.815018e-03
# [0][-14116][6351] = 3.389844e-02
# [0][-14115][6352] = 3.844029e-03
# [0][-14115][6353] = 5.366407e-03
# [0][-14115][6351] = 3.425544e-05
# [1][-14121][6351] = 3.370991e-01
# [1][-14120][6351] = 6.343077e-01
# [1][-14119][6352] = 1.641278e-01
# [1][-14119][6351] = 1.000000e+00
# [1][-14118][6352] = 1.000000e+00
# [1][-14118][6353] = 9.726192e-03
# [1][-14118][6350] = 1.532478e-02
# [1][-14118][6351] = 1.000000e+00
# [1][-14117][6352] = 1.000000e+00
# [1][-14117][6353] = 6.535208e-01
# [1][-14117][6351] = 1.000000e+00
# [1][-14116][6352] = 1.000000e+00
# [1][-14116][6353] = 7.220435e-01
# [1][-14116][6351] = 5.312908e-05
# [-1][0][0] = 0.000000e+00




# from CompactEnvelopeBuilder import PyGrid3D
# myGrid = PyGrid3D()  # Okay
# newGrid = PyGrid3D(myGrid) # Not okay


















# tfailSec = 100.
# inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
# input = open(inFileName, 'rb')
# cur_mpc = pickle.load(input)
# input.close()

# TJC.PlotDebrisFromExplodeTime(curMission, profiles, tfail=100., cutoffNAS = False)


if doMain:

    # tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles)
    # print "tProactive = {0}\n".format(tProactive)

    footprintStart = 0.
    footprintUntil = 180.

    footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)

    # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
    # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    # outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(mainFootprintFile)



if addStageReentry:
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
    firstStageFootprint = ceb.PyFootprint(outfileStr, True)

    firstStageFootprint.SmoothedOut(0)   #I believe this will simply smooth the footprints and not alter the timesteps

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
    totalFootprint = ceb.PyFootprint(mainFootprintFile, True)
    totalFootprint.MergeFootprintVectors(firstStageFootprint)
    totalFootprint.StoreFootprintAsVector(totalFootprintFile)
    totalFootprint.ExportGoogleEarth(firstStageMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)




