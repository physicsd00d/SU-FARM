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
freshWind   = False
freshDebris = False
debug       = False

doMain      = True
# addStageReentry = False


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

# from copy import deepcopy


# These parameters get injected into the final footprint name
vehicleName     = LaunchProviders.LynxMII
launchLocation  = LaunchSites.America
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
propagationParamFile = []                   # Points to thrust profile for doing propagations
precomputedParamFile = 'HTHL_Abridged.txt'  # Points to file with precomputed profile for nominal trajectory
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
curMission['debrisCatPath']           = curMission['pathToMissionFiles'] + 'DebrisCatalog/'
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
curMission['reactionTimeMinutes']       = 5     # The number of minutes that the NAS needs to safely handle a sudden debris event.
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
curMission['all_points_delta_t']      = 60.0    # Seconds, this will be the time resolution of a compact envelope
                                                #       should be GREATER THAN OR EQUAL to deltaT
curMission['numPiecesPerSample']      = 1      # The number of pieces to consider within each debris group
curMission['useAircraftDensityMap']   = False   # Do we use a uniform or the MIT density map?
curMission['debrisTimeLimitSec']      = 1*3600  # This is how long to propagate a trajectory for.  If it hasn't landed yet, then give up.

curMission['numNodes']                  = 2 # Will need to install pp to use more nodes
curMission['numNodesEnvelopes']         = 1
curMission['NASkm']                     = NASkm


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
curMission['launchAzimuth'] = 0.    #degrees, this is what it looks like on Google Earth

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

# #### ===== DEBUG =========
# # tfail = 32 has no envelope despite debris that looks like 31 and 33.
# tfailSec = 32.

# from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud#, PyFootprint
# import pickle
# import numpy as np

# deltaXY                 = curMission['deltaXY']
# deltaZ                  = curMission['deltaZ']
# h1                      = curMission['h1']
# h2                      = curMission['h2']
# debrisPickleFolder      = curMission['debrisPickleFolder']
# footprintVectorFolder   = curMission['footprintVectorFolder']
# thresh                  = curMission['thresh']
# cumulative              = curMission['cumulative']
# whichProbability        = curMission['whichProbability']

# inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
# input = open(inFileName, 'rb')
# cur_mpc = pickle.load(input)
# input.close()

# arefMeanList = cur_mpc['arefMeanList']
# numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']

# # Package them up into a PointCLoud
# # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
# curPointCloud = PyPointCloud(cur_mpc, tfailSec, curMission)

# # Place the cloud into a Grid
# curSkyGrid    = PySkyGrid(curPointCloud, deltaXY, deltaXY, deltaZ)

# # for ix in range(300):
# #     probGrid = curSkyGrid.SendGridToPython(ix)
# #     print ix
# #     if len(probGrid) > 0:
# #         print "break!"
# #         break

# # Now if I ASH without ASHing, that should just give me the unspread probabilities
# print 'ASHING'
# # h1                        = curMission['deltaXY']     # Smoothing parameters for the ASH.  Should be >= deltaXY
# # h2                        = curMission['deltaXY'] 
# h1                        = curMission['h1']     # Smoothing parameters for the ASH.  Should be >= deltaXY
# h2                        = curMission['h2'] 
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

# # tx = 299
# # print curNorm-1.

# # After uploading the density map, we have to generate the hazard probabilities
# print 'generateHazardProbabilities'
# curSkyGrid.generateHazardProbabilities(numberOfPiecesMeanList)

# STORE_IX = -666
# whichProb = 1   # probNoImpact
# hazardProbs = dict()
# for tx in range(300):
#     hazardProbs[tx] = curSkyGrid.SendProbabilitiesToPython(STORE_IX,tx,whichProb)

# # Looking at the time of failure, all the cells nearby have a nearly 1% chance of impact.  These cells should DEFINITELY get blocked off 
# #   when finding the cumulative FAA stuff.  Unless the probability of failure at this point is ridiculously low.  Must find out what
# #   pFail is getting used.  

# failProfile                 = curMission['failProfile']
# failProfileSeconds          = curMission['failProfileSeconds']
# pFail                       = curMission['pFail']
# deltaTFail                  = curMission['deltaTFail']

# ix = 0
# timelo = tfailSec
# # numGridsHere = int(np.round((timehi - timelo)/deltaTFail)) # TODO: Fix this when removing overlapping times
# # Figure out which failure times are worth propagating (i.e. they have a nonzero probability of happening)
# # timeRange = []
# # pFailThisTimestepVec = []
# # for ix in range(numGridsHere):
# sublo = timelo + ix*deltaTFail
# subhi = sublo + deltaTFail
# indices = np.where((failProfileSeconds >= sublo) & (failProfileSeconds < subhi))[0]
# pFailThisTimestep = np.sum(failProfile[indices]) * pFail

#     # if pFailThisTimestep > 0.:
#     #     timeRange.append(sublo)
#     #     pFailThisTimestepVec.append(pFailThisTimestep)

# print 'generateAllPoints_CumulativeFAA'
# curPFail = pFailThisTimestep
# EV_strike = curSkyGrid.generateAllPoints_CumulativeFAA(thresh, whichProbability, curPFail)

# haz = curSkyGrid.SendHazardPointsToPython()

# sys.exit()

# # Do the ASH
# # Must do the whole thing up-front.  On the fly only works with risk calculations at certain predetermined points.

# # 32: {0: {-18752: {7324: 0.9924537162379596,





















# sys.exit()   # Don't want to accidentally overwrite the current debris and whatnot until i squash the bugs
if debug:
    # Change a few values
    curMission['debrisCatFile']           = 'Debug.txt'
    curMission['reactionTimeMinutes']       = 5     # The number of minutes that the NAS needs to safely handle a sudden debris event.
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
    t_lo = .0
    t_hi = 180. #520.

    TJC.MonteCarlo_until_tfail(curMission, profiles, t_lo, t_hi)

# # ## Find the time until the airspace can become reactive
# minTime = 280.
# maxTime = 320.
# tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
# print "tProactive = {0}\n".format(tProactive)
# TJC.PlotNominalTrajectories(profiles, curMission, maxTime)
# sys.exit()

import numpy as np
t_lo = .0
t_hi = 179.
deltaTFail = curMission['deltaTFail']
timeVec = np.arange(t_hi*1.0,t_lo-deltaTFail,-deltaTFail)        #curTime is in seconds
for curTime in timeVec:
    print curTime
    TJC.PlotDebrisFromExplodeTime(curMission, profiles, curTime, cutoffNAS = True)
    #TJC.PlotSubEnvelopes(curMission, curTime)
sys.exit()

footprintIntervals = curMission['all_points_delta_t']
vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))
vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
mainFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '.dat'
totalFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '_stageDown.dat'

if doMain:

    # tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles)
    # print "tProactive = {0}\n".format(tProactive)

    footprintStart = 0.
    footprintUntil = 180. #520.

    footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)

    # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
    # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    # outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(mainFootprintFile)


