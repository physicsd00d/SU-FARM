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
curMission['healthMonitoringLatency'] = 0.      # Seconds

curMission['numNodes']                  = 4 # Will need to install pp to use more nodes
curMission['numNodesEnvelopes']         = 1
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

runDev = True
if runDev:
    armLength = 10.

    import numpy as np
    def getEnvelopeTimesAndFailProbs(curMission, timelo, timehi):
        deltaTFail                  = curMission['deltaTFail']
        failProfile                 = curMission['failProfile']
        failProfileSeconds          = curMission['failProfileSeconds']
        pFail                       = curMission['pFail']

        numGridsHere = int(np.round((timehi - timelo)/deltaTFail)) # TODO: Fix this when removing overlapping times

        # Figure out which failure times are worth propagating (i.e. they have a nonzero probability of happening)
        timeRange = []
        pFailThisTimestepVec = []
        for ix in range(numGridsHere):
            sublo = timelo + ix*deltaTFail
            subhi = sublo + deltaTFail
            indices = np.where((failProfileSeconds >= sublo) & (failProfileSeconds < subhi))[0]
            pFailThisTimestep = np.sum(failProfile[indices]) * pFail

            if pFailThisTimestep > 0.:
                timeRange.append(sublo)
                pFailThisTimestepVec.append(pFailThisTimestep)

        print 'times are ' + str(timeRange)
        print 'probs are ' + str(pFailThisTimestepVec)
        return timeRange, pFailThisTimestepVec



    ### Prototyping the use of non-zero health monitoring.

    # Reaction time is applied within PointCloud, so all envelopes are automatically for t <= J

    # This makes a footprint for a failure at ix over all time up to the reaction time, f + delta_R
    # genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) 
    # So that's almost the argument of the second sum, but the pointcloud will need to keep up to f + delta_R + delta_H
    # TODO: Change PointCloud to incorporate f + delta_R + delta_H

    from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud, PyGrid3D, PyFootprint

    def getProbImpacts(curMission, tfailSec):

        delta_H = int(curMission['healthMonitoringLatency']/curMission['deltaT'])    # TODO round to integer
        delta_R = int(curMission['reactionTimeSeconds']/curMission['deltaT'])  # TODO round to integer
        debrisPickleFolder      = curMission['debrisPickleFolder']

        inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
        input = open(inFileName, 'rb')
        cur_mpc = pickle.load(input)
        input.close()

        # Package them up into a PointCLoud
        curPointCloud           = PyPointCloud(cur_mpc, tfailSec, curMission)

        # Place the cloud into a Grid
        curSkyGrid              = PySkyGrid(curMission=curMission, pointCloud=curPointCloud)

        # ASH them
        # h1                        = curMission['deltaXY']     # Smoothing parameters for the ASH.  Should be >= deltaXY
        # h2                        = curMission['deltaXY'] 
        h1                        = curMission['h1']     # Smoothing parameters for the ASH.  Should be >= deltaXY
        h2                        = curMission['h2'] 
        curSkyGrid.generateASH(h1, h2)

        # Calculate all of the hazard probabilities
        curSkyGrid.generateHazardProbabilities(cur_mpc['numberOfPiecesMeanList'])

        # Finally get the probabilities for the times we want
        whichProb = curSkyGrid.getProbImpactCode()
        P_RH = curSkyGrid.GenerateSpatialProbability(whichProb, tfailSec + delta_R + delta_H, tfailSec)
        P_H = curSkyGrid.GenerateSpatialProbability(whichProb, tfailSec + delta_H, tfailSec)
        # TODO: I'm enforcing deltaTFail = 1 already, so make delta_R and delta_H just be in seconds.

        return P_RH, P_H


    # Round all the times to integers to make the dictionary lookups safer
    delta_H = int(curMission['healthMonitoringLatency']/curMission['deltaT'])    
    delta_R = int(curMission['reactionTimeSeconds']/curMission['deltaT'])  
    deltaTFail = int(curMission['deltaTFail'])
    thresh = curMission['thresh']

    probUnknown = dict()    # Holds the probability grids for the statetimes at which we don't have a VHM update yet
    probUpdated = dict()    # Holds the previous risks we were exposed to before the VHM updates told us we were safe.
    sumUpdatedProbs = PyGrid3D(); # Add in the updated probabilities as we go.

    # Flow should go like this.  Knowing delta_R and delta_H ahead of time
    # time = 0
    # SkyGrid for P_I(x, f=0 | t <= f=0 + delta_R + delta_H) to be used immediately
    # SkyGrid for P_I(x, f=0 | t <= f=0 + delta_H) to be used once the health update comes in
    footprintStart = 0.
    footprintUntil = 180.

    #TODO: MUST USE THE FAIL PROBABILITIES!!! 
    timeRange, pFailThisTimestepVec = getEnvelopeTimesAndFailProbs(curMission, footprintStart, footprintUntil)

    for tx in range(len(timeRange)):
        tfailSec = timeRange[tx]
        curPFail = pFailThisTimestepVec[tx]

        # This calculation doesn't account for the probability of failure, i.e. pFail = 1.
        P_RH, P_H = getProbImpacts(curMission, tfailSec)
        P_RH    *= curPFail                                 # So do that now
        P_H     *= curPFail

        # P_H will start to get used at curTime + delta_H + deltaTFail and will be used at every subsequent step as well.
        curTime = int(tfailSec)
        placeAtThisTime = curTime + delta_H + deltaTFail
        if placeAtThisTime <= footprintUntil:
            probUpdated[placeAtThisTime] = P_H

        # This P_RH will get used for timesteps up to curTime + delta_H
        # Need these to be their own separate objects, so use copy constructor to make them different
        tempTime = curTime
        while (tempTime <= tfailSec + delta_H):
            if tempTime <= footprintUntil:
                if tempTime in probUnknown:
                    probUnknown[tempTime] += P_RH
                else:
                    probUnknown[tempTime] = PyGrid3D(P_RH)
            tempTime += deltaTFail

        # Now take the probs at the current time, which is tfailSec, and make an envelope out of them
        if curTime in probUpdated:
            sumUpdatedProbs += probUpdated[curTime]
            del probUpdated[curTime] # remove that one from the dict since it's already been used

        # This is the total cumulative spatial probability at this time
        # P_Now = sumUpdatedProbs + probUnknown[curTime]
        P_Now = probUnknown[curTime]
        # P_Now = P_RH
        # P_Now = sumUpdatedProbs

        ## TODO: Turn off the hazard areas where the danger is completely passed
        # Don't need to block off xyz in sumUpdatedProbs if probUnknown(xyz) = 0
        # if probUnknown(xyz) == 0  ->  sumUpdatedProbs(xyz) = 0
        # if probUnknown(xyz) > 0   ->  sumUpdatedProbs(xyz) = sumUpdatedProbs(xyz)
        # sumUpdatedProbsPruned = sumUpdatedProbs.removeNoDanger(probUnknown[curTime])


        # Now that you've used it, delete it to save memory
        del probUnknown[curTime]    

        # Put this into a SkyGrid object so we can apply the threshold and make a footprint
        skyNow = PySkyGrid(curMission=curMission)
        curEV = skyNow.applyCumulativeThreshold(P_Now, thresh, np.array([int(tfailSec/deltaTFail)]))
        myFootprint = PyFootprint(skygrid=skyNow, armLength=armLength)
        myFootprint.ExportGoogleEarth(curMission['footprintVectorFolder'] + '/fpNew_' + str(tfailSec) + '.kml', yyyy, mm, dd, hour, min)

        # Store it
        outfileStr = curMission['footprintVectorFolder'] + '/fpVec_' + str(tfailSec) + '.dat'
        myFootprint.StoreFootprintAsVector(outfileStr)

        print "t = {0}, curPFail = {1}, curEV = {2}".format(tfailSec, curPFail, curEV)


        # All done here, increment time
        tfailSec += deltaTFail

    # Now that we've made them all, merge them appropriately
    # Oh boy, this is confusing.  curMission['all_points_delta_t'] is the delta_t for the final envelopes, but not for the sub-envelopes from before

    def makeFootprintForTimes(curMission, timelo, timehi):
        # For all the times [timelo, timehi], load up those vectors and merge em
        tfailSec = timelo
        while tfailSec <= timehi:
            outfileStr = curMission['footprintVectorFolder'] + '/fpVec_' + str(tfailSec) + '.dat'
            if tfailSec == timelo:
                totalFootPrint = ceb.PyFootprint(footprintFileName=outfileStr)
            else:
                totalFootPrint.MergeFootprintVectors(ceb.PyFootprint(footprintFileName=outfileStr))

            # Increment and repeat
            tfailSec += deltaTFail

        # This converts the timestep of the footprint to be the argument.  Also combines all points together and rewraps them 
        #  to produce a smooth and concise footprint.
        #totalFootPrint.SmoothedOut(newDeltaT=curMission['all_points_delta_t'], armLength=armLength)
        return totalFootPrint


        # numRange = totalFootPrint.getNumRange()
        # print "numRange = {0}".format(numRange)
        # totalFootPrint.SmoothedOut(newDeltaT=numRange * curMission['all_points_delta_t'])  # This will make footprintDelaT = numRange, and then change numRange to = 1


    for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
        timelo = footprintStart + ix*footprintIntervals
        timehi = np.min( (footprintStart + (ix+1)*footprintIntervals - deltaTFail, footprintUntil) )
        tString = "{0}-{1}".format(timelo, timehi)
        print tString

        # Make a footprint with the final timing, assembled from all the subfootprints between timelo and timehi
        curFP = makeFootprintForTimes(curMission, timelo, timehi)

        # Merge them all together
        if ix == 0:
            footprintTotal = curFP
        else:
            print '\n\nMERGE'
            footprintTotal.MergeFootprintVectors(curFP)


        # Check it out
        # curFP.ExportGoogleEarth(curMission['footprintVectorFolder'] + '/fp_' + tString + '.kml', yyyy, mm, dd, hour, min)


    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    sys.exit()

# If still debugging
## TODO: Check the curPFails for each time
## TODO: Check the returned values of risk for each sub envelope


debugSingleTime = False
if debugSingleTime:
     ### ===== DEBUG =========
    tfailSec = 70.

    from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud#, PyFootprint
    import pickle
    import numpy as np

    deltaXY                 = curMission['deltaXY']
    deltaZ                  = curMission['deltaZ']
    h1                      = curMission['h1']
    h2                      = curMission['h2']
    debrisPickleFolder      = curMission['debrisPickleFolder']
    footprintVectorFolder   = curMission['footprintVectorFolder']
    thresh                  = curMission['thresh']
    cumulative              = curMission['cumulative']
    whichProbability        = curMission['whichProbability']

    inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
    input = open(inFileName, 'rb')
    cur_mpc = pickle.load(input)
    input.close()

    arefMeanList = cur_mpc['arefMeanList']
    numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']

    # Package them up into a PointCLoud
    # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
    curPointCloud = PyPointCloud(cur_mpc, tfailSec, curMission)

    # Place the cloud into a Grid
    curSkyGrid    = PySkyGrid(curMission=curMission, pointCloud=curPointCloud)

    # # Now if I ASH without ASHing, that should just give me the unspread probabilities
    print 'ASHING'
    h1                        = curMission['deltaXY']     # Smoothing parameters for the ASH.  Should be >= deltaXY
    h2                        = curMission['deltaXY'] 
    # h1                        = curMission['h1']     # Smoothing parameters for the ASH.  Should be >= deltaXY
    # h2                        = curMission['h2'] 
    curSkyGrid.generateASH(h1, h2)

    def checkNorm(ash):
        curNorm = 0.
        for curZ in ash:
            for curX in ash[curZ]:
                for curY in ash[curZ][curX]:
                    curNorm += ash[curZ][curX][curY]
        return curNorm

    # Okay, now I can look through the histograms any way I want
    # curID = 10  # highest beta
    curID = 2   # most pieces.  This must SURELY generate a hazard area.  Very light, mostly hangs in air.
    hist = dict()
    ash = dict()
    whichProb = 0   # Impact
    for tx in range(300):
        hist[tx] = curSkyGrid.SendHistogramToPython(curID,tx)
        ash[tx] = curSkyGrid.SendProbabilitiesToPython(curID,tx, 0)

        # if len(hist[tx]) > 0:
        # print "{0}: {1} --> {2}".format(tx, hist[tx], ash[tx])
        print "{0}: {1} --> {2}".format(tx, hist[tx], 1-checkNorm(ash[tx]))





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




