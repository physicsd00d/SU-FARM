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
# # Points to the python scripts needed from Francisco
# friscoFiles = '../../../Prop3Dof/FriscoDebris/pythonFiles/'
# # Points to the binaries for propagating trajectories
# debrisPropPATH = '../../../Prop3Dof/FriscoDebris/'
# # Points to the files that I've written
# tjcFiles = '../../'
#
# import os
# import sys
# sys.path.append(friscoFiles)
# sys.path.append(debrisPropPATH)
# sys.path.append(tjcFiles)

# Build directory
buildDir = "../../../../build/"
tjcFiles = '../../' # This is where the TJC.py file resides.

import os
import sys
# sys.path.append(friscoFiles)
# sys.path.append(debrisPropPATH)
sys.path.append(os.path.abspath(buildDir))
sys.path.append(os.path.abspath(tjcFiles))


'''
Import the modules necessary for the script.
Note: Must Come After sys.path update
Note: This block is probably the same across all main scripts
'''
# import debrisPythonWrapper as dpw
# import getPropTraj as traj
# import AtmosProfile as AP
import CompactEnvelopeBuilder as ceb
# import numpy as np
import TJC
import datetime as dt

import LaunchSites
import LaunchProviders



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
# from ReadMission import readInput

# Initialize the mission
propagationParamFile = []                   # Points to thrust profile for doing propagations
                                            #   If you use this option, then you need to worry about dtval!
precomputedParamFile = 'HTHL_Abridged.txt'  # Points to file with precomputed profile for nominal trajectory
pathToMissionFiles = './'                   # Kind of a holdover from a previous file structure

# Planet info
omegaE = 7.2921158494529352e-05             # rad/s
planetModel = 0                             # 0 means spherical, 1 means elliptical

curMission = dict(propagationParamFile = propagationParamFile, precomputedParamFile = precomputedParamFile,
                pathToMissionFiles = pathToMissionFiles, omegaE = omegaE, planetModel = planetModel)
curMission = TJC.InitializeMission(curMission)


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
curMission['debrisCatFile']           = 'LynxDebrisCatalog.txt'
# curMission['debrisCatFile']           = 'columbiaWithBlast.txt'
# curMission['atmosphere']              = friscoFiles + 'AtmoProfiles/special.txt'
# curMission['atmospherePickle'] = '../AtmoProfiles/FrontRange.pkl'
curMission['atmospherePickle'] = '../AtmoProfiles/SpaceportAmerica.pkl'


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
ExportDate = dt.datetime(year=yyyy, month=mm, day=dd, hour=hour, minute=min, second=sec)

curMission['ExportDate'] = [yyyy, mm, dd, hour, min]
curMission['ExportDateDT'] = ExportDate

'''
#
# #
# ======================= Begin Computations ============================ #
# #
#
'''

########### Fold that into existing infrastructure

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
    # t_hi = 200.
    # t_lo = 305.0

    t_hi = 520.
    # t_hi = 240.

    TJC.MonteCarlo_until_tfail(curMission, profiles, t_lo, t_hi)

minTime = 280.
maxTime = 320.
tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
print "tProactive = {0}\n".format(tProactive)
sys.exit()


if not debug:
    footprintStart = 0.
    # footprintUntil = 110.
    footprintUntil = 520.
    # footprintUntil = 240.
    footprintIntervals = 60.

    footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)
    vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))

    # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
    # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(outfileStr)





else:
#     print 'SPECIAL DEBUGGING'
#
#     # Expand this function with copy paste so we can dissect it a little
#     # EV_strike1, zeroToOneFootprint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
#     # Know this from normal debugging
#     # times are [115.0]
#     # probs are [0.0009803921550000001]
#     tfailSec = 145.0
#     curPFail = 0.0009803921550000001
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
#     EV_strike = curSkyGrid.generateAllPoints(numberOfPiecesMeanList, arefMeanList, thresh, curPFail)
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
#     # // The FAA uses a rather coarse grid, so convert to their coarse grid
#     # //      These are the parameters of the coarsened grid
#     newDeltaXY   = 3.5      #//[km]
#     newDeltaZ    = 20.      #//[km]  This is higher than NASkm, but I need the values to nest, so hopefully this is fine
#     spatialProb = curSkyGrid.GetSpatialProbabilty_Coarse(newDeltaXY, newDeltaZ)
#
#
#     # Transform into something we can plot
#     xyProbVector = []
#     # deltaXY = curMission['deltaXY']
#     deltaXY = newDeltaXY
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
#     plt.savefig("GeneratedFiles/GraphLynx{0}.png".format(deltaXY))
#
#
# # x, y = np.random.rand(2, 100) * 4
# # hist, xedges, yedges = np.histogram2d(x, y, bins=4)
# #
# # elements = (len(xedges) - 1) * (len(yedges) - 1)
# # xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)
# #
# # xpos = xpos.flatten()
# # ypos = ypos.flatten()
# # zpos = np.zeros(elements)
# # dx = 0.5 * np.ones_like(zpos)
# # dy = dx.copy()
# # dz = hist.flatten()
# #
# # ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
# #
# # plt.show()
#
#
#     sys.exit()













    # print 'NORMAL DEBUGGING'
    # timelo = 115.
    # # timelo = 50.
    # timehi = 115.
    # TJC.PlotDebrisFromExplodeTime(curMission, profiles, timehi*1.0)
    # totalNumTimeSteps = -1; #Don't use this anymore
    #
    # numEventsSimulated = profiles['numTrajSamples'] * profiles['numWindSamples']     # do this properly in the future, save into MPC
    # numDebrisPerIXSimulated = numEventsSimulated * curMission['numPiecesPerSample']   # This will give us the normalization for the probability.
    # numDebrisPerIXSimulated = 0
    #
    # print 'footprint'
    # # Let's make the updating footprint for zeroToOne minutes
    # EV_strike1, zeroToOneFootprint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
    #
    # # ExportDate = dt.datetime(year=yyyy, month=mm, day=dd, hour=hour, minute=min, second=sec)
    #
    # startTimeMinutes = ExportDate.hour*60. + ExportDate.minute
    # offsetTimeMinutes = 0
    # tstepMinutes = int(curMission['all_points_delta_t']/60.)
    # # TODO: Change the tstep units to seconds
    # zeroToOneFootprint.MakeFacetFiles(curMission['facetFolder'], startTimeMinutes, offsetTimeMinutes, tstepMinutes)
    # print 'done'
    # sys.exit()



    # print 'DEBUGGING NEW WAY OF GENERATING ENVELOPES'
    # # Set the footprint timestep to 1 second because that's the easiest to conceptualize
    # curMission['all_points_delta_t']      = 1.0    # Seconds, this will be the time resolution of a compact envelope
    # curMission['thresh']                    = 1e-7
    # curMission['deltaZ']                    = NASkm/4.   #km
    # curMission['deltaTFail']              = 5.0     # Seconds, this is how often we explode the rocket

    # tfailSec = 100.0

    # EVstrike, Fprint = TJC.makeFootprintFromTimes(curMission, tfailSec, tfailSec)
    # Let's start by looking at the first 10 seconds for the very first fail time
    # Load the precomputed footprint


    # ''' Prototype Flash'''
    footprintStart = 0.
    # # footprintUntil = 50.
    footprintUntil = 510.
    footprintIntervals = 60.
    # footprintTotal = []

    # TJC.GenerateEnvelopes_Flash(curMission, footprintStart, footprintUntil, footprintIntervals)

    #
    # for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
    #     timelo = footprintStart + ix*footprintIntervals
    #     timehi = np.min( (footprintStart + (ix+1)*footprintIntervals, footprintUntil) )
    #
    #     print 'TIMES: From {0} to {1}'.format(timelo, timehi)
    #     EVstrike, curFootPrint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
    #     print 'EV =  ' + str(EVstrike)
    #
    #     # Now take that footprint and...
    #     # Smooth it out to a single timestep
    #     numRange = curFootPrint.getNumRange()
    #     curFootPrint.SmoothedOut(numRange)  # This will make footprintDelaT = numRange, and then change numRange to = 1
    #
    #     numRange = curFootPrint.getNumRange()
    #     FPDeltaT = curFootPrint.getDeltaT()
    #
    #     # Resize the deltaT to be only the length of the interval
    #     #   So if we're making an envelope at each second, then the footprint should be chopped at 1 second
    #     #   If we're combining times, like every 5 seconds or every minute, then it should be 5s or 60s
    #     curFootPrint.ChopTimeAt(footprintIntervals)
    #
    #     # Translate the footprint forward to tfailSec
    #     if timelo > 0:
    #         curFootPrint.SlideFootprintBySeconds(timelo)
    #
    #     # Merge it with the others
    #     if ix == 0:
    #         footprintTotal = curFootPrint
    #     else:
    #         print '\n\nMERGE'
    #         footprintTotal.MergeFootprintVectors(curFootPrint)
    #
    #
    #     # Print to GE
    #     debugFolder = 'GeneratedFiles/Sandbox/'
    #     vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, timelo, timehi)
    #     curFootPrint.ExportGoogleEarth(debugFolder + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    #
    #     # # Fprint.SmoothedOut(footprintIntervals)
    #     #
    #     # # Fprint.SmoothedOut()
    #     # curFootPrint.ExportGoogleEarth('GeneratedFiles/PythonGE_' + str(timelo) + 'To'
    #     #                                       + str(timehi) + 'FootprintSMOOTH.kml', yyyy, mm, dd, hour, min)
    #
    #
    #
    #
    # # Just to be safe(?), set the params we need in order to translate / rotate
    # footprintTotal.SetAzimuthDeg(curMission['launchAzimuth'])
    # footprintTotal.SetLaunchLatDeg(curMission['launchLat'])
    # footprintTotal.SetLaunchLonDeg(curMission['launchLon'])
    #
    # # Print to GE
    # debugFolder = 'GeneratedFiles/Sandbox/'
    # vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, footprintStart, footprintUntil)
    # # footprintTotal.ExportGoogleEarth(debugFolder + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    #
    # footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    #
    # outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    # footprintTotal.StoreFootprintAsVector(outfileStr)

    sys.exit()







    # Load the footprint
    curFPVecFile = '{0}/fpVec_{1}.dat'.format(curMission['footprintVectorFolder'], tfailSec)
    curFootPrint = ceb.PyFootprint(curFPVecFile, True)

    # Smooth it out to a single timestep
    numRange = curFootPrint.getNumRange()
    curFootPrint.SmoothedOut(numRange)  # This will make footprintDelaT = numRange, and then change numRange to = 1

    numRange = curFootPrint.getNumRange()
    FPDeltaT = curFootPrint.getDeltaT()

    # Resize the deltaT to be only 1 second long
    curFootPrint.ChopTimeAt(1.)

    # Translate the footprint forward to tfailSec
    curFootPrint.SlideFootprintBySeconds(tfailSec)

    # Merge it with the others

    # Print to GE
    debugFolder = 'GeneratedFiles/Sandbox/'
    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, tfailSec, numRange)
    curFootPrint.ExportGoogleEarth(debugFolder + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)

    # Exit
    sys.exit()




    # Keep only the first 10 seconds
    fpLengthSec = -1
    # fpLengthSec = curFootPrint.ChopTimeAt(curMission['reactionTimeMinutes'] * 60.)

    # Print to GE
    debugFolder = 'GeneratedFiles/Sandbox/'
    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, tfailSec, fpLengthSec)
    curFootPrint.ExportGoogleEarth(debugFolder + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)

    # Open up the debris
    input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
    cur_mpc = pickle.load(input)
    input.close()

    from data2GE import convertTJC

    debrisOutFile = '{0}debris_{1}.kml'.format(debugFolder, tfailSec)
    flatArray = cur_mpc['flatPointArray']
    numTimeSteps = cur_mpc['numTimeSteps']
    numRuns = len(numTimeSteps)
    convertTJC(debrisOutFile, flatArray, numTimeSteps, numRuns, cutoffNAS = False, maxTimeSteps = 1e10)


    myTraj = ceb.PyTrajectory()
    myTraj.loadDebrisTrajectory(cur_mpc, tfailSec, curMission)
    myTraj.ExportGoogleEarth(debrisOutFile, ExportDate)























# # outfileStr = GeneratedFilesFolder + 'LynxMII_SSA.dat'
# outfileStr = 'OtherPythonFiles/FootprintLibrary/LynxMII_SSA.dat'
# thisFP = ceb.PyFootprint(outfileStr, True)
# thisFP.SetAzimuthDeg(curMission['launchAzimuth'])
# thisFP.StoreFootprintAsVector(outfileStr)






# ## Find the time until the airspace can become reactive
# maxTime = 510
# tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, maxTime)
# print 'tProactive = ' + str(tProactive)
# sys.exit()

# print 'DEBUG'
# timelo = 0.
# #timelo = 50.
# timehi = timelo
# TJC.PlotDebrisFromExplodeTime(curMission, profiles, timelo*1.0)
# totalNumTimeSteps = -1; #Don't use this anymore
#
# numEventsSimulated = profiles['numTrajSamples'] * profiles['numWindSamples']     # do this properly in the future, save into MPC
# numDebrisPerIXSimulated = numEventsSimulated * curMission['numPiecesPerSample']   # This will give us the normalization for the probability.
# numDebrisPerIXSimulated = 0
#
# print 'footprint'
# # Let's make the updating footprint for zeroToOne minutes
# EV_strike1, zeroToOneFootprint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
# print 'done'
# sys.exit()
