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
# Points to the python scripts needed from Francisco
friscoFiles = '../../../Prop3Dof/FriscoDebris/pythonFiles/'
# Points to the binaries for propagating trajectories
debrisPropPATH = '../../../Prop3Dof/FriscoDebris/'
# Points to the files that I've written
tjcFiles = '../../'

import os
import sys
sys.path.append(friscoFiles)
sys.path.append(debrisPropPATH)
sys.path.append(tjcFiles)


'''
Import the modules necessary for the script.
Note: Must Come After sys.path update
Note: This block is probably the same across all main scripts
'''
import debrisPythonWrapper as dpw
import getPropTraj as traj
import AtmosProfile as AP
import CompactEnvelopeBuilder as ceb
import numpy as np
import TJC
import datetime as dt

import LaunchSites
import LaunchProviders

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
# from ReadMission import readInput

# Initialize the mission
propagationParamFile = 'nominalParam_new.txt'   # Points to thrust profile for doing propagations
precomputedParamFile = []                       # Points to file with precomputed profile for nominal trajectory
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
# curMission['debrisCatFile']           = 'testFile.txt'
curMission['debrisCatFile']           = 'Halcon9_1stNEW.txt'
# curMission['debrisCatFile']           = 'columbiaWithBlast.txt'
curMission['atmospherePickle'] = '../AtmoProfiles/Cape.pkl'


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

# Precompute some Atmosphere and Trajectory profiles
freshWind   = False
freshDebris = False
debug       = False

doMain      = False
addStageReentry = True

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
    t_hi = 180.

    TJC.MonteCarlo_until_tfail(curMission, profiles, t_lo, t_hi)

# # ## Find the time until the airspace can become reactive
# minTime = 120.
# maxTime = 180.
# tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
# print "tProactive = {0}\n".format(tProactive)


# TJC.PlotNominalTrajectories(profiles, curMission, 180.)
# sys.exit()

footprintIntervals = curMission['all_points_delta_t']
vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))
vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
mainFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '.dat'
totalFootprintFile = curMission['footprintLibrary'] + vehicleFileName + '_stageDown.dat'

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

    # sys.exit()


if addStageReentry:
    '''# Prototype handling the first stage reentry'''

    # Specify the time of the staging
    tStage = 175

    curMission['debrisPickleFolder']      = curMission['GeneratedFilesFolder']  + 'debrisPickleFolder'

    firstStageMission = deepcopy(curMission)
    firstStageMission['debrisPickleFolder'] = curMission['GeneratedFilesFolder']  + 'firstStagePickleFolder'
    firstStageMission['debrisCatFile'] = 'firstStage.txt'
    firstStageMission['reactionTimeMinutes'] = -1
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














#
# if not debug:
#
#     # tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles)
#     # print "tProactive = {0}\n".format(tProactive)
#
#     footprintStart = 0.
#     footprintUntil = 180.
#     footprintIntervals = 60.
#
#     footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)
#     vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))
#
#     # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
#     # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))
#
#     vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)
#
#     footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
#     outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
#     footprintTotal.StoreFootprintAsVector(outfileStr)
#
#     sys.exit()
#












# else:
#     print 'SPECIAL DEBUGGING'
#
#     # Expand this function with copy paste so we can dissect it a little
#     # EV_strike1, zeroToOneFootprint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
#     # Know this from normal debugging
#     # times are [115.0]
#     # probs are [0.0009803921550000001]
#     tfailSec = 45.0
#     curPFail = 1./180.
#
#     #
#     # Then it calls genFootprint(curMission, tfailSec, curPFail), expanded here
#     #
#
#     ExportDateDT            = curMission['ExportDateDT']
#     # reactionTimeMinutes     = curMission['reactionTimeMinutes']       # Where did this go?
#     deltaXY                 = curMission['deltaXY']
#     deltaZ                  = curMission['deltaZ']
#     h1                      = curMission['h1']
#     h2                      = curMission['h2']
#     debrisPickleFolder      = curMission['debrisPickleFolder']
#     footprintVectorFolder   = curMission['footprintVectorFolder']
#     thresh                  = curMission['thresh']
#     useAircraftDensityMap   = curMission['useAircraftDensityMap']       # Do we use a uniform or the MIT density map?
#     cumulative              = curMission['cumulative']
#     whichProbability        = curMission['whichProbability']
#
#     from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud, PyFootprint
#
#     # Open up the debris
#     input = open(debrisPickleFolder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#     cur_mpc = pickle.load(input)
#     input.close()
#
#     arefMeanList = cur_mpc['arefMeanList']
#     numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']
#
#     # This is [total number of pieces simulated within this mpc] / [number of debris categories in this mpc]
#     # TODO: If all_points_delta_t != debrisDeltaT, then we'll be double-counting here.
#     # numDebrisPerIXSimulated = cur_mpc['numPieces']/len(numberOfPiecesMeanList)
#
#     # Package them up into a PointCLoud
#     # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
#     curPointCloud = PyPointCloud(cur_mpc, tfailSec, curMission)
#
#     # Place the cloud into a Grid
#     curSkyGrid    = PySkyGrid(curPointCloud, deltaXY, deltaXY, deltaZ)
#
#     # Do the ASH
#     # Must do the whole thing up-front.  On the fly only works with risk calculations at certain predetermined points.
#     curSkyGrid.generateASH(h1, h2)
#
#     # print 'tfailsec = ' + str(tfailSec)
#     # return -1, './GeneratedFiles/footprintVectorFolder/fpVec_20.0.dat'
#
#     # Get the lat/lons of the filled cells
#     latlonArray = curSkyGrid.createEmptyAircraftDensityMap()
#
#     if useAircraftDensityMap:
#         # With those lat/lons, find the probAircraft for each cell
#         from AircraftDensityMap import AircraftDensityMap as ADM
#         density = ADM()
#         densityArray = density.getDensity(latlonArray, ExportDateDT)
#         print 'THIS IS PROBABLY A DENSITY AND NOT A PROBABILITY.  FIX THIS!!!'
#         sys.exit()
#
#         # Send that information back into C++
#         curSkyGrid.populateAircraftDensityMap(densityArray, len(densityArray))
#     else:
#         fourNM2 = 13.72                     #// 4 (n.m.)^2 * (1.852 km/nm)^2 = 13.72 km^2
#         aircraftDensity = 1./fourNM2        #// [prob/km^2] Paul Wilde's assumed aircraft density (1 every 4nm^2)
#         cellArea = deltaXY*deltaXY
#         # print '\n\n\nDEBUG: SETTING PROB AIRPLANE TO ONE\n\n\n'
#         probOfAirplaneInCell = aircraftDensity * cellArea;
#
#         import numpy as np
#         # Send that information back into C++
#         curSkyGrid.populateAircraftDensityMap(np.array([probOfAirplaneInCell]), -1)
#
#     # EV_strike = curSkyGrid.generateAllPoints(numberOfPiecesMeanList, arefMeanList, thresh, curPFail)
#     # After uploading the density map, we have to generate the hazard probabilities
#     curSkyGrid.generateHazardProbabilities(numberOfPiecesMeanList, curPFail)
#
#     # Now apply the current definition of 'cumulative' for the chosen type of probability
#     if cumulative == 'FAA':
#         # # This will require coarsening the grid
#         # newDeltaXY   = 3.5      #//[km]
#         # newDeltaZ    = 20.      #//[km]  This is higher than NASkm, but I need the values to nest, so hopefully this is fine
#         #
#         # print 'SHIT GUYS!
#         # raise RuntimeError
#         newDeltaXY   = -1      #//[km]
#         newDeltaZ    = -1
#
#         # Inside this function, the grid will be permanently coarsened.  I don't think it's a good idea to try
#         #   to alter the curMission here because if it's shared memory than that could potentially mess up other
#         #   calculations that are happening in parallel threads.  Wait until whole thing is over and then do it.
#         EV_strike = curSkyGrid.generateAllPoints_CumulativeFAA(thresh, whichProbability, newDeltaXY, newDeltaZ)
#
#     elif cumulative == 'TJC':
#         EV_strike = curSkyGrid.generateAllPoints_CumulativeTJC(thresh, whichProbability)
#
#
#
#
#
#
#     #
#     # NOW WE CAN DEBUG!!!  Retrieve the spatial probability map, then plot it
#     #
#
#     # PROB_IMPACT      = 1001
#     # PROB_CASUALTY    = 1002
#     # PROB_CATASTROPHE = 1003
#     # curSkyGrid.GenerateSpatialProbability(PROB_IMPACT)
#
#     # This is a dict of dicts of dicts
#     # spatialProb = curSkyGrid.GetSpatialProbabilty()
#
#     # # // The FAA uses a rather coarse grid, so convert to their coarse grid
#     # # //      These are the parameters of the coarsened grid
#     # newDeltaXY   = 3.5      #//[km]
#     # newDeltaZ    = 20.      #//[km]  This is higher than NASkm, but I need the values to nest, so hopefully this is fine
#     # spatialProb = curSkyGrid.GetSpatialProbabilty_Coarse(newDeltaXY, newDeltaZ)
#
#     spatialProb = curSkyGrid.GetSpatialProbabilty()
#
#
#     # Transform into something we can plot
#     xyProbVector = []
#     # deltaXY = curMission['deltaXY']
#     # deltaXY = newDeltaXY
#
#     zVals = spatialProb.keys()
#     curZ = zVals[0]
#     for curX in spatialProb[curZ].keys():
#         for curY in spatialProb[curZ][curX].keys():
#             curProb = spatialProb[curZ][curX][curY]
#
#             # Have to transform the x and y's from indices into actual coordinates
#             #   Actually, doing this here will lead to bugs later on when comparing float values.
#             #   Leave as indices until the very end.
#             # xyProbVector.append([curX*deltaXY, curY*deltaXY, curProb])
#             xyProbVector.append([curX, curY, curProb])  #can't mix types in a vector, so curX and curY will be set to float
#
#     # Convert to numpy for slicing, then slice out the information
#     xyProbVector    = np.array(xyProbVector)
#     xValues         = xyProbVector[:,0].astype(int)    # These are actualy integers and you know it.
#     yValues         = xyProbVector[:,1].astype(int)
#     probValues      = xyProbVector[:,2]
#
#     # BUG.  So there appears to be some rounding error coming from somewhere.  I don't understand where
#     # Anyways, need to clean the position values up, round them all to the same precision as deltaXY
#
#     # Find the range of X and create a vector of the EXPANDED x-values of the histogram
#     minXVal     = np.min(xValues)
#     maxXVal     = np.max(xValues)
#     xEdges      = range(minXVal, maxXVal+1)
#     # numXSteps   = int((maxXVal - minXVal)/deltaXY) + 1
#     # xEdges      = np.linspace(minXVal, maxXVal, numXSteps)
#
#     # Find the range of Y and create a vector of the EXPANDED y-values of the histogram
#     minYVal     = np.min(yValues)
#     maxYVal     = np.max(yValues)
#     yEdges      = range(minYVal, maxYVal+1)
#     # numYSteps   = int((maxYVal - minYVal)) + 1
#     # yEdges      = np.linspace(minYVal, maxYVal, numYSteps)
#
#     # Inflate the sparse map into a full meshgrid
#     xpos, ypos = np.meshgrid(xEdges, yEdges)
#
#     # bar3d takes 1D arrays, so flatten everything
#     # These are the locations (presumable lower-left) from which each point is referenced
#     xpos = xpos.flatten()
#     ypos = ypos.flatten()
#     zpos = np.zeros(len(xpos))
#
#     # These are the measurements that define the size of each bar
#     dx = deltaXY * np.ones_like(zpos)
#     dy = dx.copy()
#     dz = np.zeros_like(dy)
#
#     # Determine the density threshold
#     # weightedThresh = curMission['thresh'] * (deltaXY**2) / (4* fourNM2)
#     weightedThresh = curMission['thresh']# * (deltaXY**2) / (4* fourNM2)
#     print "weightedThresh = " + str(weightedThresh)
#
#     # Keep track of the area that gets blocked off
#     affectedArea = 0.
#
#     # Now we have to fit the sparse data into this expanded mesh
#     for row in xyProbVector:
#         [curX, curY, curProb] = row
#         goodX = (xpos == curX)                      # Find the good X indices
#         goodY = (ypos == curY)                  # Find the good Y indices
#         goodIndex = np.all([goodX, goodY], 0)   # Find the one index that satisfies them both
#
#         if sum(goodIndex) != 1:
#             print '\n\nERROR'
#             print sum(goodIndex)
#             print row
#         else:
#             # Apply the threshold
#             if curProb > weightedThresh:
#                 dz[goodIndex] = weightedThresh
#                 affectedArea += 1
#             else:
#                 dz[goodIndex] = curProb
#
#             # # Don't apply the threshold
#             # dz[goodIndex] = curProb
#
#     print 'affectedArea = ' + str(affectedArea * deltaXY**2)
#     # Okay, NOW convert the xy values from indices into their true values
#     xpos = xpos*deltaXY
#     ypos = ypos*deltaXY
#
#     import matplotlib as mpl
#     mpl.use('Agg')  # Allows plot generation on server without X-windows
#
#     from mpl_toolkits.mplot3d import Axes3D
#     import matplotlib.pyplot as plt
#     from matplotlib.colors import LogNorm
#
#     # Make the colors
#     norm = LogNorm(1e-14, 1e-6)
#     colors = plt.cm.jet(norm(dz))
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#
#     # ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
#     ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average', color=colors)
#     # ax.view_init(elev=90.,azim=-90.)
#     # plt.show()
#     # This only works for wx backend
#     # See this thread: http://stackoverflow.com/questions/12439588/how-to-maximize-a-plt-show-window-using-python
#     # mng = plt.get_current_fig_manager()
#     # mng.frame.Maximize(True)
#     # mng.full_screen_toggle()
#     plt.savefig("GeneratedFiles/Graph{0}{1}.png".format(vehicleName,deltaXY))