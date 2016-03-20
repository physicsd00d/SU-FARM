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
debug                   = False

import os
import sys

curFilePath = os.path.dirname(os.path.abspath(__file__)) + "/"
rootDir =   os.path.abspath(curFilePath + "../../../../") + "/"
outputDir = rootDir + "outputs/" # Where to store results, gitignored
tempDir =   rootDir + "temp/"   # temp files here, gitignored


'''
Import the modules necessary for the script.
Note: Must Come After sys.path update
Note: This block is probably the same across all main scripts
'''
# import debrisPythonWrapper as dpw
# import getPropTraj as traj
# import AtmosProfile as AP
import CompactEnvelopeBuilder as ceb
import numpy as np
from Simulation import TJC
import datetime as dt

import matplotlib
matplotlib.use('Agg')  # Allows plot generation on server without X-windows
import matplotlib.pyplot as plt

from Simulation import LaunchSites
from Simulation import LaunchProviders

# These parameters get injected into the final footprint name
vehicleName     = LaunchProviders.Reentry
launchLocation  = LaunchSites.OK    # NOTE: even though it says 'launch', in this context it really means 'landing'
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
curMission['debrisCatPath']           = curMission['pathToMissionFiles'] + 'DebrisCatalog/'
curMission['debrisCatFile']           = 'testFileDistributed.txt'
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
curMission['reactionTimeMinutes']       = 5     # The number of minutes that the NAS needs to safely handle a sudden debris event.
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
curMission['numPiecesPerSample']      = [10]      # The number of pieces to consider within each debris group
                                                #       IF EMPTY, that means use the actual number for each debris group
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

secondsFromMidnightUTC = hour*3600 + min*60 + sec


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

profiles = []
if (freshWind):
    # Should really move all the important mission stuff into this if-statement and wrap it up into the montecarlo dictionary
    
    numTrajSamples = 1      # If you change this to anything other than 1, it might break.  Look at numDebrisPerIXSimulated to start.
    numWindSamples = 30

    # I only need to generate wind profiles here, since i'm not going to worry about multiple nominal trajectories yet
    # Could / should probably anticipate doing it though andjust replicate the single trajectory here to conform with the existing infrastrcture

    atmStorage, stateVecStorage, thetagStorage, tfailStorage\
        = TJC.GenerateWindTrajProfilesDirectional(curMission, numTrajSamples, numWindSamples, angleLow, angleHi, windMagCoeff)

    profiles = dict(atmStorage = atmStorage, stateVecStorage = stateVecStorage, thetagStorage = thetagStorage,
                    tfailStorage = tfailStorage, numTrajSamples = numTrajSamples, numWindSamples = numWindSamples)

    import pickle
    output = open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl', 'wb')
    pickle.dump(profiles,output)
    output.close()

    # sys.exit()
else:
    import pickle
    profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))


if freshDebris:
    # Generate the debris
    coeffIX = []
    lowerBreakLimit = [] # By setting these to 0 and 1, we'll explode at just the lower time
    upperBreakLimit = [] # IF YOU CHANGE THIS!!!  Then you'll need to fix the risk calculations later that assume zero
                         # If you set them both to [], then will explode at all times
    TJC.MonteCarlo_Distributed_Reentry_Wrapper_CAIB(curMission, coeffIX, curMission['numPiecesPerSample'],
                                                    lowerBreakLimit, upperBreakLimit, profiles)
    # TJC.MonteCarlo_Distributed_Reentry_Wrapper_CAIB(curMission, coeffIX, curMission['numPiecesPerSample'], profiles)


    # # We need to first make sure that our initial state vector is far enough back in the trajectory to be useful
    # # How long from initialStateVec until debris hits the NAS?
    # numWindIX = len(profiles['atmStorage'])
    # print 'numWindIX = ' + str(numWindIX)
    #
    # minTimeMinutes = 1e9
    # for windIX in range(numWindIX):
    #
    #     outFileName = '{0}/mpc_{1}_{2}.pkl'.format(curMission['debrisPickleFolder'], windIX, int(tfail))
    #
    #     input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'rb')
    #     cur_mpc = pickle.load(input)
    #     input.close()
    #
    #     minTimeMinutes = np.min((minTimeMinutes, TJC.FindShortestTimeAboveNAS(curMission, cur_mpc)))
    #
    #     # '''DEBUG'''
    #     # # Now ship this off to a google earth routine
    #     # import data2GE
    #     # GEfile = '../../GeneratedFiles/GE_Debris.kml'
    #     # data2GE.convertTJC(GEfile, cur_mpc['flatPointArray'], cur_mpc['numTimeSteps'], len(cur_mpc['numTimeSteps']))
    #
    # print minTimeMinutes



# windIX = 6
# input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'rb')
# cur_mpc = pickle.load(input)
# print TJC.FindShortestTimeAboveNAS(curMission, cur_mpc)
# input.close()

# sys.exit()





## Find the time until the airspace can become reactive
#tProactive = TJC.FindStateTimeForProactiveArchitecture(mission1, profiles)



if not debug:
    curMission['numWindProfiles'] = len(profiles['atmStorage'])

    footprintStart = profiles['tfailStorage'][0][0][0]
    footprintUntil = profiles['tfailStorage'][0][0][-1]
    # footprintStart = 240.
    # footprintUntil = 300.
    footprintIntervals = 60.

    print 'footprintStart = {0}'.format(footprintStart)
    print 'footprintUntil = {0}'.format(footprintUntil)

    footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)
    vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))

    # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
    # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(outfileStr)

else:
    print 'DEBUGGING'


























if plotColumbiaGround:
    # ==== Now Make A Plot Of The Debris ====
    # Use this to plot
    finalLat = []
    finalLon = []

    numWindIX = len(profiles['atmStorage'])
    print 'numWindIX = ' + str(numWindIX)

    lowerBreakLimit = profiles['tfailStorage'][0][0][0] #Assuming first trajectory is representative
    upperBreakLimit = profiles['tfailStorage'][0][0][-1]
    deltaTFail = curMission['deltaTFail']
    # Won't include very last time step which may well be negative altitude and thus undesireable
    timeVec = np.arange(lowerBreakLimit,upperBreakLimit,deltaTFail)        #curTime is in seconds

    for windIX in range(numWindIX):
        for tfail in timeVec:
            inFileName = '{0}/mpc_{1}_{2}.pkl'.format(curMission['debrisPickleFolder'], windIX, int(tfail))
            # input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'rb')
            input = open(inFileName, 'rb')
            cur_mpc = pickle.load(input)
            input.close()
            sizeFlatPointArray = cur_mpc['sizeFlatPointArray']
            numTimeSteps = cur_mpc['numTimeSteps']
            flatPointArray = cur_mpc['flatPointArray']

            for ix in range(len(numTimeSteps)):
                curDebris = flatPointArray[ (6*np.sum(numTimeSteps[:ix])) : (6*np.sum(numTimeSteps[:(ix+1)])) ].reshape((numTimeSteps[ix],6))
                finalLat.append(curDebris[-1][0])
                finalLon.append(curDebris[-1][1])

    # for plotting purposes
    im = plt.imread('Sandbox/debrisField.png')
    latShift = 0.086     # Increasing this number moves the test point downward
    lonShift = -0.01      # Increasing this number moves the test point leftward

    leftLon     = -98.2617 + lonShift
    rightLon    = -92.1090 + lonShift
    bottomLat   = 29.6854 + latShift -0.04
    upperLat    = 34.0093 + latShift
    corners = [leftLon,rightLon,bottomLat,upperLat]
    # corners = [-98.2617,-92.1090,29.6854,34.0093]
    implot = plt.imshow(im,extent=corners)

    plt.scatter(finalLon, finalLat,marker='o', color='red', s=.1)
    plt.savefig('../../GeneratedFiles/ResultsMain')
    # sys.exit()

    # ==== Here's snippets that relate to the more general distributed breakup ====
    # lowerBreakLimit = 15.
    # upperBreakLimit = 135.
    # lowerBreakLimit = -1
    # upperBreakLimit = -1
    # TJC.MonteCarlo_Distributed_Reentry_Wrapper(mission1, coeffIX, mission1['numPiecesPerSample'], profiles)
    # cur_mpc = TJC.MonteCarlo_Distributed_Reentry(mission1, coeffIX, mission1['numPiecesPerSample'],
    #                                              profiles, lowerBreakLimit, upperBreakLimit, windIX)
    # cur_mpc = TJC.MonteCarlo_Distributed_Reentry(mission1, [], 5, profiles)

if calcIndividualHazard:
    # ============ First going to load up the aircraft tracks ============

    # Read in the aircraft information
    input = open('HighRiskOutput.txt', 'r')
    acid2typeMap = dict()
    curTrackTime = 0

    ft2km       = 0.0003048
    knots2km_s  = 0.000514444444
    # Going to leave the top layer as a dict / map because acids are not necessarily going to start at zero or appear in order
    # Be careful also to realize that will a python dict will let you mix types in a layer, std::map will not
    aircraftRecord = dict()
    for line in input:
        key = line.split()

        if len(key) > 0:
            # There's something here
            if (len(key) == 1) and (int(key[0]) < 24*3600):
                # This is a track time
                curTrackTime = int(key[0])
            elif (len(key) == 6):
                # This is an aircraft track
                acid = int(key[0])
                acType = key[1]
                curLat = float(key[2])
                curLon = float(key[3])
                curLevel = float(key[4]) * 100. * ft2km
                curSpeed = float(key[5]) * knots2km_s

                # Attempt a conversion from geodetic to spherical
                # [curLat, curLon, _] = TJC.geodetic2spherical(curLat,curLon,curLevel*1e3)


                # Make sure python knows the structure if this is the first time seeing the acid
                if not aircraftRecord.has_key(acid):
                    aircraftRecord[acid] = [],acType    # It's a tuple, which will get translated to a std::pair

                aircraftRecord[acid][0].append([float(curTrackTime), curLat, curLon, curLevel, curSpeed])
    input.close()


    # ============ Now load the debris, grid it, and ASH it ============

    # Get the debris
    import pickle
    profiles = pickle.load(open(pathToMissionFiles + 'localProfiles.pkl','rb'))
    numDebrisPickles = profiles['numWindSamples']


    ### Eventually this will all get wrapped into a function
    import CompactEnvelopeBuilder as ceb

    curSkyGrid = []
    # Load each of the pickle files into one large skygrid
    for windIX in range(numDebrisPickles):
    # for windIX in range(2):
        input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'rb')
        cur_mpc = pickle.load(input)
        input.close()

        # Package them up into a PointCLoud
        timeOfInitialFailure = 0
        curPointCloud = ceb.PyPointCloud(cur_mpc, timeOfInitialFailure, curMission)     # Applies all_points_delta_t

        # Loading the pointcloud throws out points past the reaction time, still need to chop???
        # newNumSteps = curPointCloud.ChopAfterSeconds(int(tfailSec + 60*reactionTimeMinutes))
        #print '(' + str(tfailSec) + ',' + str(newNumSteps) + ')'

        # Place the cloud into a Grid
        if windIX == 0:
            curSkyGrid    = ceb.PySkyGrid(curPointCloud, curMission['deltaXY'], curMission['deltaXY'], curMission['deltaZ'])
        else:
            curSkyGrid.IncorporatePointCloudIntoGrid(curPointCloud)

    # Making the assumption that all profiles use the exactly same debris catalog
    arefMeanList = cur_mpc['arefMeanList']
    numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']

    # Upload the aircraft data (locations, areas, type)
    curSkyGrid.UploadAircraft(aircraftRecord)

    from AircraftAreaModel import AircraftAreaModel as AAM
    curSkyGrid.UploadAircraftPropertiesMap(AAM)

    #Assuming all debris categories had same number of pieces simulated
    numUniqueIDs = len(np.unique(cur_mpc['debrisID']))
    numTotalPieces = cur_mpc['numPieces']

    # Basically, this is a count of how many individual pieces of debris were simulated, within each
    #  ballistic coefficient group, over the all_points timestep.  So if you simulate a single piece of debris with
    #  a timestep of 1 sec, but then grid the sky in 60 second increments, it's like you simulated that piece 60 times
    #  over that larger timestep.  Put another way, you've double-counted that same piece 60 times when histogramming.
    #
    # The above is only true if the debris timestep is 1.0.  Really, you're double-counting by all_points_delta_t/debrisDeltaT
    # TODO: Unfortunately, that time step is not an option.  Attend to this soon
    # TODO: Make these changes in TJC.genFootprint as well
    #
    # numDebrisPerIXSimulated = numDebrisPickles * numTotalPieces * (mission1['all_points_delta_t'])/numUniqueIDs
    # print '\n\n\nWHEN YOU ARE LESS TIRED, CHECK THAT MULTIPLYING BY DELTA T IS THE RIGHT THING TO DO HERE!!!!!\n\n\n'

    # # ASH the probabilities around
    # curSkyGrid.generateASH(mission1['h1'], mission1['h2'], numDebrisPerIXSimulated)
    #
    # numDebrisTimeSteps = curSkyGrid.getNumRange()
    # print "numDebrisTimeSteps = {0}".format(numDebrisTimeSteps)
    #
    # # Do risk calculation
    # riskToAC = curSkyGrid.CalculateRiskToIndividualAircraft(numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC)


    riskToAC = curSkyGrid.CalculateRiskToIndividualAircraft_OnTheFly(numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC,
                                                   curMission['h1'], curMission['h2'])


    for acid in riskToAC:
        print '{0} = {1}'.format(aircraftRecord[acid][1], riskToAC[acid])

    print '\nCorrecting for difference in timesteps'
    for acid in riskToAC:
        print '{0} = {1}'.format(aircraftRecord[acid][1], riskToAC[acid]/curMission['all_points_delta_t'])





    '''
    This is for making animations.  It's long and convoluted, sorry.
    If you don't want one, you can stop the program here!
    '''
    makeAnimation = True
    if not makeAnimation:
        sys.exit()

    from scipy.interpolate import griddata
    # import matplotlib
    # matplotlib.use('Agg')  # Allows plot generation on server without X-windows
    # import matplotlib.pyplot as plt

    import matplotlib.animation as animation
    import types # used in duck punching

    np.set_printoptions(linewidth=120)

    from matplotlib.colors import LogNorm
    import matplotlib.gridspec as gridspec

    # for plotting purposes
    fig1 = plt.figure()


    latShift = 0.086     # Increasing this number moves the test point downward
    lonShift = -0.01      # Increasing this number moves the test point leftward
#
# # for plotting purposes
# im = plt.imread('debrisField.png')
# corners = [-98.2617 + lonShift, -92.1090 + lonShift, 29.6854 + latShift -0.04, 34.0093 + latShift]

    debrisField = plt.imread('Sandbox/debrisField.png')
    leftLon     = -98.2617 + lonShift
    rightLon    = -92.1090 + lonShift
    bottomLat   = 29.6854 + latShift -0.04
    upperLat    = 34.0093 + latShift
    corners = [leftLon,rightLon,bottomLat,upperLat]
    # corners = [-98.2617,-92.1090,29.6854,34.0093]
    # implot = ax1.imshow(debrisField,extent=corners)


    # TIME LOOP WILL START HERE
    altIX2plot  = 0
    lowAlt      = altIX2plot * curMission['deltaZ']
    highAlt     = (altIX2plot+1) * curMission['deltaZ']

    # ASH the probabilities around
    print "ASHING"
    curSkyGrid.generateASH(curMission['h1'], curMission['h2'])

    # numDebrisTimeSteps will each be all_points_delta_t long
    numDebrisTimeSteps = curSkyGrid.getNumRange()
    print "numDebrisTimeSteps = {0}".format(numDebrisTimeSteps)

    figCounter = 1
    images = []
    for tx in range(numDebrisTimeSteps):
        # Get the individual probability of impact grid
        probGrid = curSkyGrid.SendGridToPython(tx)
        if len(probGrid) == 0:
            continue    # Debris hasn't made it to the NAS yet or there are no planes in the area of the debris

        # The indices that correspond to the lowest alt level
        altitudeLevels = np.unique(probGrid[:,0])
        if len(altitudeLevels) < (altIX2plot+1):
            continue    # no probabilities to plot here.  move on
        altIX = (probGrid[:,0] == altitudeLevels[altIX2plot])
        altGrid = probGrid[altIX,:]

        # Find the limits of the current grid in order to interpolate over them
        minLat = np.min(altGrid[:,2])
        maxLat = np.max(altGrid[:,2])
        minLon = np.min(altGrid[:,1])
        maxLon = np.max(altGrid[:,1])

        # This is the mesh to interpolate over based on above limits
        lonVec = np.linspace(minLon,maxLon,100)
        latVec = np.linspace(minLat,maxLat,100)

        # probInterp = griddata((altGrid[:,1], altGrid[:,2]), altGrid[:,3], (lonVec[None,:], latVec[:,None]), method='cubic')
        probInterp = griddata((altGrid[:,1], altGrid[:,2]), altGrid[:,3], (lonVec[None,:], latVec[:,None]), method='linear')

        print 'tx = {0} with a timestep of {1}'.format(tx, curMission['all_points_delta_t'])


        # plt.clf()
        # plt.imshow(debrisField,extent=corners)   # Hopefully this will wipe the old
        # single_im = plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
        # plt.text(0,1,"tx = {0:04d}".format(tx))

        # # duck punching so that contour works
        # def setvisible(self,vis):
        #    for c in self.collections: c.set_visible(vis)
        # single_im.set_visible = types.MethodType(setvisible,single_im,None)
        # single_im.axes = plt.gca()
        # images.append([single_im])

        plt.clf()
        gs = gridspec.GridSpec(1, 2,width_ratios=[15,1])
        lowerBoundExponent = -4
        upperBoundExponent = 0

        # do the 2nd subplot, the pseudo colorbar, first
        ax2 = plt.subplot(gs[1])
        # np.logspace gives you logarithmically spaced levels -
        # this, however, is not what you want in your colorbar
        #
        # you want equally spaced labels for each exponential group:
        #
        levls = np.array([])
        for lix in range(lowerBoundExponent, upperBoundExponent):
            levls = np.concatenate((levls, np.linspace(10**lix,10**(lix+1),10)))

        # levls = np.linspace(1,10,10)
        # levls = np.concatenate((levls[:-1],np.linspace(10,100,10)))
        # levls = np.concatenate((levls[:-1],np.linspace(100,1000,10)))
        # levls = np.concatenate((levls[:-1],np.linspace(1000,10000,10)))

        #
        # simple x,y setup for a contourf plot to serve as colorbar
        #
        XC = [np.zeros(len(levls)), np.ones(len(levls))]
        YC = [levls, levls]
        CM = ax2.contourf(XC,YC,YC, levels=levls, norm = LogNorm())
        # log y-scale
        ax2.set_yscale('log')
        # y-labels on the right
        ax2.yaxis.tick_right()
        # no x-ticks
        ax2.set_xticks([])



        ax1 = plt.subplot(gs[0])
        ax1.imshow(debrisField,extent=corners)   # Hopefully this will wipe the old
        # plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
        # cbarTicks = np.linspace(0, 1, 11, endpoint=True)
        # cbarTicks = np.log10(np.logspace(1e-5, 1, endpoint=True))
        cbarTicks = np.logspace(-4,0)


        ax1.contourf(lonVec,latVec,probInterp, cbarTicks, cmap=plt.cm.jet, norm=LogNorm(), vmin=0., vmax=1.)
        # ax1.colorbar(ticks=cbarTicks)
        # plt.colorbar(ticks=cbarTicks, spacing='proportional')

        # print 'cbarTicks = ' + str(cbarTicks)
        plt.title("tx = {0:04d}, dt = {1} sec, altRange = [{2}, {3}]km".format(tx, curMission['all_points_delta_t'], lowAlt, highAlt))

        # Find the AC tracks at this time
        for key in aircraftRecord:
            curAC = np.array(aircraftRecord[key][0])
            # goodIX = np.where((curAC[:,0] - secondsFromMidnightUTC) == 1.*(tx))
            # print 'goodIX = ' + str(goodIX)

            acType = aircraftRecord[key][1]

            # goodIX = (curAC[:,0] - secondsFromMidnightUTC) == 1.*(tx)*mission1['all_points_delta_t']
            timesHere       = range(tx*int(curMission['all_points_delta_t']), (tx+1)*int(curMission['all_points_delta_t']) )
            goodIX          = np.in1d(curAC[:,0].ravel() - secondsFromMidnightUTC, timesHere).reshape(curAC[:,0].shape)
            totalACLat      = curAC[goodIX,1]
            totalACLon      = curAC[goodIX,2]
            totalACLevel    = curAC[goodIX,3]

            # print 'times here = ' + str(timesHere)
            # print 'len = ' + str(len(totalACLat))
            # Is there anybody here?
            # if len(totalACLat) > 0:
            for pt in range(len(totalACLat)):
                curLat = totalACLat[pt]
                curLon = totalACLon[pt]
                curAlt = totalACLevel[pt]
                # print '{0}, {1}'.format(curLon, curLat)
                if ((upperLat > curLat) & (curLat > bottomLat) & (leftLon < curLon) & (curLon < rightLon) & (lowAlt <= curAlt) & ( curAlt < highAlt)):
                #     print 'HERE!!!!'
                    ax1.scatter(curLon ,curLat,marker='o',s=5, color='red')

                # If last pt
                if pt == (len(totalACLat)-1):
                    ax1.annotate(acType, xy=(totalACLon[pt],totalACLat[pt]))

        plt.savefig("_tmp{0:04d}.png".format(figCounter))
        # images.append( (plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k'),) )
        figCounter += 1

        # plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
        # plt.contourf(lonVec,latVec,probInterp,15,cmap=plt.cm.jet)

    im_ani = animation.ArtistAnimation(fig1, images, interval=1, repeat_delay=1, blit=True)

    # im_ani = animation.Animation(fig1, images, blit=True)
    # im_ani.save('GeneratedFiles/basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    im_ani.save('doesntmatter.mp4', clear_temp=False)
    from subprocess import call
    call(["ffmpeg -y -r 5 -i _tmp%04d.png -b:v 1800k ../../GeneratedFiles/Contours.mp4"], shell=True)
    call(["mv _tmp*.png MoviePngs/"], shell=True)



    # plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
    # plt.savefig('GeneratedFiles/PlotGridProb')

    print 'exiting'
    sys.exit()











































# #numGrids = int(np.round((testUpper - testLower)/deltaTFail) + 1)
# #probWeights = 1./numGrids
#
# #thresh = 1e-7
# numEventsSimulated = profiles['numTrajSamples'] * profiles['numWindSamples']     # do this properly in the future, save into MPC
# numDebrisPerIXSimulated = numEventsSimulated * mission1['numPiecesPerSample']   # This will give us the normalization for the probability.
#
# ## Let's make the updating footprint for zeroToOne minutes
# #timelo = 5.
# #timehi = 55.
# #EV_strike1, zeroToOneFootprint = TJC.makeFootprintFromTimes(mission1, deltaTFail, thresh, totalNumTimeSteps, numDebrisPerIXSimulated, timelo, timehi)
# #print 'returned ' + str(EV_strike1)
# ##zeroToOneFootprint.ChopTimeAt(1)
# ##zeroToOneFootprint.ShiftFootprintByMinutes(2)   #This is additional minutes, for for a 3 minute footprint, want to shift by 3-1=2
# #zeroToOneFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_zeroToOneFootprint.kml', yyyy, mm, dd, hour, min)
# #zeroToOneFootprint.SmoothedOut()
# #zeroToOneFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_zeroToOneFootprintSMOOTH.kml', yyyy, mm, dd, hour, min)
# #
# ## Let's make the updating footprint for oneToTwo minutes
# #timelo = 60.
# #timehi = 115.
# #EV_strike2, oneToTwoFootprint = TJC.makeFootprintFromTimes(mission1, deltaTFail, thresh, totalNumTimeSteps, numDebrisPerIXSimulated, timelo, timehi)
# #print 'returned ' + str(EV_strike2)
# ##oneToTwoFootprint.ChopTimeAt(1)
# ##oneToTwoFootprint.ShiftFootprintByMinutes(2)
# #
# #oneToTwoFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_oneToTwoFootprint.kml', yyyy, mm, dd, hour, min)
# #oneToTwoFootprint.SmoothedOut()
# #oneToTwoFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_oneToTwoFootprintSMOOTH.kml', yyyy, mm, dd, hour, min)
# #
# #sys.exit()
# #print 'made it to last one'
#
# # Let's make the updating footprint for twoToThree minutes
# timelo = 120.
# timehi = 170.
# EV_strike3, twoToThreeFootprint = TJC.makeFootprintFromTimes(mission1, deltaTFail, thresh, totalNumTimeSteps, numDebrisPerIXSimulated, timelo, timehi)
# print 'returned ' + str(EV_strike3)
#
# #twoToThreeFootprint.ChopTimeAt(1)
# #twoToThreeFootprint.ShiftFootprintByMinutes(2)
#
# twoToThreeFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_twoToThreeFootprint.kml', yyyy, mm, dd, hour, min)
# twoToThreeFootprint.SmoothedOut()
# twoToThreeFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_twoToThreeFootprintSMOOTH.kml', yyyy, mm, dd, hour, min)
# sys.exit()
#
#
# # Make the facet stuff
# zeroToOneFootprint.MergeFootprintVectors(oneToTwoFootprint)
# zeroToOneFootprint.MergeFootprintVectors(twoToThreeFootprint)
# zeroToOneFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_totalFootprint.kml', yyyy, mm, dd, hour, min)
#
#
# #exit()
# #
# #folderName = 'GeneratedFiles/testSUA'
# #startTimeMinutes = hour*60 + min
# #offsetTimeMinutes = 0
# #tstepMinutes = mission1['all_points_delta_t']
# #zeroToOneFootprint.MakeFacetFiles(folderName, startTimeMinutes, offsetTimeMinutes, tstepMinutes)
#
#
#
# ##
# ##











































#for ix in range(numGrids):
#    tfailSec = testLower + ix*deltaTFail
#    print 'tfailSec = ' + str(tfailSec)
#    input = open(folder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#    cur_mpc = pickle.load(input)
#    
#    input.close()
#    
#    curPointCloud = ceb.PyPointCloud(cur_mpc, tfailSec)
#    curSkyGrid    = ceb.PySkyGrid(curPointCloud, mission1['deltaXY'], mission1['deltaXY'], 5)
#    curSkyGrid.generateAllPointsFromSimpleHistogram(thresh, totalNumTimeSteps, numEventsSimulated)
#
#    myFootprint = ceb.PyFootprint(curSkyGrid)
#    if (ix == 0):
#        totalFootprint = ceb.PyFootprint(curSkyGrid)
#    else:
#        totalFootprint.MergeFootprintVectors(myFootprint)
#        
#    #    myFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_' + str(tfailSec) + '.kml', yyyy, mm, dd, hour, min)
#
##totalFootprint.SmoothedOut()
#totalFootprint.ProjectAllPointsDown()
#totalFootprint.ChopTimeAt(tProactive/60.)
#totalFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_totalFootprint.kml', yyyy, mm, dd, hour, min)



















#print "\n\n  ~~~~~~~ MARCHING FORWARD ~~~~~~~~~~~ \n\n"
#
## Initialize the skygrid
#tfailSec = testLower
#
#input = open(folder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#cur_mpc = pickle.load(input)
#input.close()
#
## PointCloud expects the deltaT to be in minutes
#curPointCloud = ceb.PyPointCloud(cur_mpc, tfailSec)
#curSkyGrid    = ceb.PySkyGrid(curPointCloud, mission1['deltaXY'], mission1['deltaXY'], 5)
#
#tfailSec += deltaTFail
#
#while (tfailSec <= testUpper):
#    input = open(folder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#    cur_mpc = pickle.load(input)
#    input.close()
#
#    curPointCloud = ceb.PyPointCloud(cur_mpc, tfailSec)
#    #curSkyGrid.EncorporateDebrisIntoGrid( cur_mpc['flatPointArray'], cur_mpc['numTimeSteps'],  mission1, cur_mpc['UTC'], tfailSec)
#    curSkyGrid.IncorporatePointCloudIntoGrid( curPointCloud)
#    print 'just completed ' + str(tfailSec)
#
#    tfailSec += deltaTFail
#
#curSkyGrid.generateAllPointsFromGrid()
#myFootprint = ceb.PyFootprint(curSkyGrid)
#myFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_forward.kml', yyyy, mm, dd, hour, min)



# == Loop
# Explode at each timestep
#tfail = tProactive
#TJC.MonteCarlo_until_tfail(mission1, profiles, tfail, folder)

## BONUS POINTS ##
# Bin the points
# ADVANCED: Do ASH or KDE


#input = open(folder + '/mpc_' + str(55.0) + '.pkl', 'rb')
#cur_mpc = pickle.load(input)
#input.close()
#
#
## Write them to file for posterity
## == Loop Ends
#curPointCloud = ceb.PyPointCloud(cur_mpc)
#curSkyGrid    = ceb.PySkyGrid(curPointCloud, .5, .5, .5)
#curSkyGrid.ConvertToProbability()
#
##input = open(folder + '/mpc_' + str(45.0) + '.pkl', 'rb')
##cur_mpc = pickle.load(input)
##input.close()
##
##curPointCloud = ceb.PyPointCloud(cur_mpc)
##curSkyGrid.inco
#
#
##curSkyGrid.ExportToGoogleEarth(mission1['GoogleEarthFile'], yyyy, mm, dd, hour, min)
#
## Add em together in some way that makes sense




#deltaTFail = 5.0
#
## Initialize the skygrid
#tfailSec = 120.0
#
#input = open(folder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#cur_mpc = pickle.load(input)
#input.close()
#
## PointCloud expects the deltaT to be in minutes
#curPointCloud = ceb.PyPointCloud(cur_mpc, tfailSec)
#curSkyGrid    = ceb.PySkyGrid(curPointCloud, mission1['deltaXY'], mission1['deltaXY'], 5)
#
#tfailSec -= deltaTFail
#
#while (tfailSec > 0):
#    input = open(folder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#    cur_mpc = pickle.load(input)
#    input.close()
#
#    curPointCloud = ceb.PyPointCloud(cur_mpc, tfailSec)
#    #curSkyGrid.EncorporateDebrisIntoGrid( cur_mpc['flatPointArray'], cur_mpc['numTimeSteps'],  mission1, cur_mpc['UTC'], tfailSec)
#    curSkyGrid.EncorporatePointCloudIntoGrid( curPointCloud)
#    print 'just completed ' + str(tfailSec)
#
#    tfailSec -= deltaTFail
#
#curSkyGrid.generateAllPointsFromGrid(mission1['deltaXY'])
#myFootprint = ceb.PyFootprint(curSkyGrid)
#myFootprint.ExportGoogleEarth('GeneratedFiles/PythonGE_backward.kml', yyyy, mm, dd, hour, min)







## Verify the trajectory
#tfail = profiles['tfailStorage']
#stateVecs = profiles['stateVecStorage']
##stateIX = np.where(tfail[0][0] == 120)[0][0]
#Vx = stateVecs[0][0][3]
#Vy = stateVecs[0][0][4]
#Vz = stateVecs[0][0][5]
#vel = np.sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
#
#
#stateIX = np.where(tfail[0][0] == 55)[0][0]
#print 'At tfail = ' + str(tfail[0][0][stateIX]) + '~~~~~~~~~~~~~'
#print 'alt  = ' + str(latlonaltStorage[3*stateIX+2]/1000) + ' km'
#print 'vel  = ' + str(vel[stateIX]/1000) + ' km/s'
#print 'down = ' + str(TJC.HaversineDistance(latlonaltStorage[:2], latlonaltStorage[3*stateIX:(3*stateIX + 2)])) + ' km'
#print '~~~~~~~~~~~~~~\n\n\n'
#
#stateIX = np.where(tfail[0][0] == 115)[0][0]
#print 'At tfail = ' + str(tfail[0][0][stateIX]) + '~~~~~~~~~~~~~'
#print 'alt  = ' + str(latlonaltStorage[3*stateIX+2]/1000) + ' km'
#print 'vel  = ' + str(vel[stateIX]/1000) + ' km/s'
#print 'down = ' + str(TJC.HaversineDistance(latlonaltStorage[:2], latlonaltStorage[3*stateIX:(3*stateIX + 2)])) + ' km'
#print '~~~~~~~~~~~~~~\n\n\n'
#
#stateIX = np.where(tfail[0][0] == 145)[0][0]
#print 'At tfail = ' + str(tfail[0][0][stateIX]) + '~~~~~~~~~~~~~'
#print 'alt  = ' + str(latlonaltStorage[3*stateIX+2]/1000) + ' km'
#print 'vel  = ' + str(vel[stateIX]/1000) + ' km/s'
#print 'down = ' + str(TJC.HaversineDistance(latlonaltStorage[:2], latlonaltStorage[3*stateIX:(3*stateIX + 2)])) + ' km'
#print '~~~~~~~~~~~~~~\n\n\n'





