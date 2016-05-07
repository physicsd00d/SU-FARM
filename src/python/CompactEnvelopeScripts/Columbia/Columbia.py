'''
#
This file calculates the risk to aircraft from the Columbia accident
#
'''

'''
System imports that let python know about the file structure of the code
Note: This block should be the same across all main scripts

If you want to use valgrind to debug, compile with -g for debug symbols then call
valgrind --tool=memcheck --suppressions=valgrind-python.py python -E -tt falcon9.py
Be sure to remove -g from compilation when done otherwise code will be slooooow
'''
freshMain   			= False  # State Vector
freshWind               = True  # Why uncertain wind for this case? B/c uncertainty in direction is manually tweaked.
freshDebris             = True
debug                   = False

plotColumbiaGround      = False
calcIndividualHazard    = False
makeAnimation           = False     # Turn this off if using a small all_pts_delta_t
                                    # ASH grid must be pretty coarse for this to work, but why?

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

# TODO: Add a columbia reentry thing
# These parameters get injected into the final footprint name
vehicleName     = LaunchProviders.Reentry
launchLocation  = LaunchSites.OK    # NOTE: even though it says 'launch', in this context it really means 'landing'
vehicleNotes    = ct.vehicleNotes + '_Columbia'

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
from TrajectoryInfo import initVec

# Initialize the mission
propagationParamFile = []                   # Points to thrust profile for doing propagations
precomputedParamFile = 'ColumbiaMainPieceTrajectory.txt'  # Points to file with precomputed profile for nominal trajectory
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
curMission['debrisCatPath']     = debrisPath + 'Columbia/'
# curMission['debrisCatFile']           = 'testFileDistributed.txt'
# curMission['debrisCatFile']           = 'debugDistributed.txt'
curMission['debrisCatFile']           = 'debugColumbiaMarch16.txt'
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

curMission['deltaXY']                   = 10.    #km
curMission['deltaZ']                    = NASkm/1.   #km
curMission['h1']                        = 20.    # Smoothing parameters for the ASH.  Should be >= deltaXY
curMission['h2']                        = 20.

# Parameters for the safety architecture of the NAS
curMission['reactionTimeSeconds']       = -5*60.     # The number of seconds that the NAS needs to safely handle a sudden debris event.
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
curMission['deltaTFail']              = 1.0     # Seconds, this is how often we explode the rocket
curMission['all_points_delta_t']      = 60.0     # Seconds, this will be the time resolution of a compact envelope
                                                #       should be GREATER THAN OR EQUAL to deltaT
                                                #       For reentry, appears to control the deltaT of the movies made
curMission['numPiecesPerSample']      = [10]      # The number of pieces to consider within each debris group
                                                #       IF EMPTY, that means use the actual number for each debris group
curMission['useAircraftDensityMap']   = False   # Do we use a uniform or the MIT density map?
curMission['debrisTimeLimitSec']      = 1*3600  # This is how long to propagate a trajectory for.  If it hasn't landed yet, then give up.

curMission['numNodes']                  = 2
curMission['numNodesEnvelopes']         = -1 	# This should not get used.
curMission['NASkm']                     = NASkm



# '''
# Do I still need these? If they're uncommented, they're asked for in monte carlo function
# '''

# # Updated, but not sure i need these (this is the first point from the nominal trajecotry file, read this in eventually)
curMission['launchLat'] = launchLat
curMission['launchLon'] = launchLon
curMission['launchAlt'] = 0.   #km

# # Not updated, also not sure i need these
curMission['initialUTC'] = 156.84861111111113 	# (i think) This number could be anything, as long as it's consistent
												# Asked for in the monte carlo functions
# # curMission['launchAzimuth'] = 169.    #degrees, this is the heading angle of the SSA runway.  Measured with Google Earth
curMission['launchAzimuth'] = 0.    #degrees, this is what it looks like on Google Earth

'''
OUTPUT options
'''
# Export files as GoogleEarth or FACET
curMission['exportGE']    = True
curMission['exportFACET'] = False

# This date gets used for GE.  I believe the FACET files are date agnostic, but i think the time of day might get set here
# This specific time of 13:59:30 is the time that the CAIB used for the beginning of the progressive breakup.
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
# Wind angle measured from East, positive is counterclock
curMission['noWind'] = True
# angleLow    = 0.
# angleHi     = 0.
angleHi         = -5.
angleLow        = -35. #-90
windMagCoeff    = 1.0   # Scale the wind profile

''' Some comments about the initial state vector'''
# For Columbia, this is the moment of initial breakup (when we have the statevector for)
# 8:59:37 EST -- TJC: the time when the shuttle lost hydraulic pressure and would have spun out of control
# (13:59 UTC)

# Entry Interface (EI) over the Pacific Ocean at 8:44:09 a.m. EST  (13:44:09 UTC)
# EI + 555 is first reports of people seeing debris being shed (13:53:09 UTC)

# 13:59:30 GMT is the time used for simulation in CAIB as being the beginning of the 2-minute breakup
# Francisco's state vector starts 15 seconds before that, so the time is 13:59:15

# 01
# Correcting for difference in timesteps
# angleHi         = -10.
# angleLow        = -30. #-90
# windMagCoeff    = 1.0   # Scale the wind profile
# E120 = 7.56701529558e-06
# MD82 = 2.76489574755e-05
# B733 = 1.49465958892e-05
# B735 = 5.54765289564e-05
# MD90 = 0.000121642723396
# B738 = 1.78835542652e-05
# CRJ2 = 0.000113516714622
# CRJ2 = 7.72987187524e-05
# B733 = 6.31568551159e-05

# 02
# Correcting for difference in timesteps
# angleHi         = -5.
# angleLow        = -35. #-90
# windMagCoeff    = 1.0   # Scale the wind profile
# E120 = 7.23672775471e-06
# MD82 = 2.6999740812e-05
# B733 = 1.46709281366e-05
# B735 = 5.31696191339e-05
# MD90 = 0.00012085102269
# B738 = 2.49021304948e-05
# CRJ2 = 0.00011157188999
# CRJ2 = 7.9430048788e-05
# B733 = 6.23311270627e-05

## This is with the new interpolated aircraft tracks.  Still need to check on the risk calculation.  1 alt level, 0.5 vs 3, nowind
# Correcting for difference in timesteps
# E120 = 0.000176470285798
# MD82 = 0.00187881155333
# B733 = 0.00104933947812
# B735 = 0.0013383864785
# MD90 = 0.00269132952572
# B738 = 0.0
# CRJ2 = 0.002097949861
# CRJ2 = 0.00144988348205
# B733 = 0.000894275342716
# Number of Pieces     = 78849
# Total Mass of Pieces = 26403.8302996

# ======================= Begin Computations ============================ #
# Precompute some Atmosphere and Trajectory profiles

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

    # It doesn't seem like te tfailStorage here would be useful.  Those are the statetimes of the vehicle's trajectory.
    # What I want to do with Columbia is explode it at a single time, possibly the first time, and then have a
    # progressive breakup.  So the number of failure timesteps would be equal to the time steps in the debris catalog.
    import pickle
    output = open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl', 'wb')
    pickle.dump(profiles,output)
    output.close()
else:
    import pickle
    profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))


if freshDebris:
    print "Make sure this works with     TJC.MonteCarloDebris(curMission, profiles, t_lo, t_hi)"
    sys.exit()  
    # Generate the debris
    coeffIX = []
    lowerBreakLimit = 0 # By setting these to 0 and 1, we'll explode at just the lower time
    upperBreakLimit = 1 # IF YOU CHANGE THIS!!!  Then you'll need to fix the risk calculations later that assume zero
    TJC.MonteCarlo_Distributed_Reentry_Wrapper_CAIB(curMission, coeffIX, curMission['numPiecesPerSample'],
                                                    lowerBreakLimit, upperBreakLimit, profiles)




if plotColumbiaGround:
    # ==== Now Make A Plot Of The Debris ====
    # Use this to plot
    finalLat = []
    finalLon = []

    numWindIX = len(profiles['atmStorage'])
    print 'numWindIX = ' + str(numWindIX)

    for windIX in range(numWindIX):

        '''This makes the assumption that we're ONLY using the zeroth timestep'''
        inputFile = '{0}/mpc_{1}_{2}.pkl'.format(curMission['debrisPickleFolder'], windIX, int(0))
        input = open(inputFile, 'rb')
        # input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'rb')
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
    im = plt.imread('debrisField.png')
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
    plt.savefig(curMission['GeneratedFilesFolder'] + "ResultsMain")
    # sys.exit()



if calcIndividualHazard:
    # ============ First going to load up the aircraft tracks ============

    # Read in the aircraft information
    # input = open('HighRiskOutput.txt', 'r')
    input = open('HighRiskETMSInterp.txt', 'r')
    aircraftTrackDeltaTSec = 1
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

    # Want to know the latest time that AC are in the air
    maxAircraftTrackTime = 0
    for (curTrack, acName) in aircraftRecord.itervalues():
        maxAircraftTrackTime = max(curTrack[-1][0],maxAircraftTrackTime)
        # print acName, (curMaxTime - secondsFromMidnightUTC)/curMission['all_points_delta_t']



    # ============ Now load the debris, grid it, and ASH it ============

    # Get the debris
    import pickle
    profiles = pickle.load(open(curMission['GeneratedFilesFolder'] + 'localProfiles.pkl','rb'))
    numDebrisPickles = profiles['numWindSamples']


    ### Eventually this will all get wrapped into a function
    # import CompactEnvelopeBuilder as ceb

    curSkyGrid = []
    # Load each of the pickle files into one large skygrid
    for windIX in range(numDebrisPickles):
    # for windIX in range(2):
        '''This makes the assumption that we're ONLY using the zeroth timestep'''
        inputFile = '{0}/mpc_{1}_{2}.pkl'.format(curMission['debrisPickleFolder'], windIX, int(0))
        input = open(inputFile, 'rb')
        # input = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'rb')
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
    arefMeanList            = cur_mpc['arefMeanList']
    numberOfPiecesMeanList  = cur_mpc['numberOfPiecesMeanList']
    debrisMass              = cur_mpc['debrisMass']
    # Upload the aircraft data (locations, areas, type)
    curSkyGrid.UploadAircraft(aircraftRecord, aircraftTrackDeltaTSec)

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

    riskToAC = curSkyGrid.CalculateRiskToIndividualAircraft_OnTheFly(numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC,
                                                   curMission['h1'], curMission['h2'])


    for acid in riskToAC:
        print '{0} = {1}'.format(aircraftRecord[acid][1], riskToAC[acid])

    print '\nCorrecting for difference in timesteps'
    for acid in riskToAC:
        print '{0} = {1}'.format(aircraftRecord[acid][1], riskToAC[acid]/curMission['all_points_delta_t'])

    print 'Number of Pieces     = {0}'.format(np.sum(numberOfPiecesMeanList))
    print 'Total Mass of Pieces = {0}'.format(np.sum(debrisMass))



















    # ============ End of Risk Calculation =====================

    plotRisk = False

    '''
    This is for making animations.  It's long and convoluted, sorry.
    If you don't want one, you can stop the program here!
    '''
    if not makeAnimation:
        sys.exit()

    from scipy.interpolate import griddata
    # import matplotlib
    # matplotlib.use('Agg')  # Allows plot generation on server without X-windows
    # import matplotlib.pyplot as plt

    # import matplotlib.animation as animation
    # import types # used in duck punching

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

    debrisField = plt.imread('debrisField.png')
    leftLon     = -98.2617 + lonShift
    rightLon    = -92.1090 + lonShift
    bottomLat   = 29.6854 + latShift -0.04
    upperLat    = 34.0093 + latShift
    corners = [leftLon,rightLon,bottomLat,upperLat]
    # corners = [-98.2617,-92.1090,29.6854,34.0093]
    # implot = ax1.imshow(debrisField,extent=corners)

    # For simplicity, ensure that deltaT = 1
    if (curMission['deltaT'] != 1.):
        print "ERROR: When making columbia movie, require debris delta_t = 1 second, not {0} seconds".format(curMission['deltaT'])
        sys.exit()

    # TIME LOOP WILL START HERE
    altIX2plot  = 0
    lowAlt      = altIX2plot * curMission['deltaZ']
    highAlt     = (altIX2plot+1) * curMission['deltaZ']

    if plotRisk:
        # ASH the probabilities around
        print "ASHING plotRisk"
        curSkyGrid.generateASH(curMission['h1'], curMission['h2'])
        print "DONE ASHING plotRisk"

    numDebrisTimeSteps = curSkyGrid.getNumRange()       # These are in curMission['deltaT']
    print "numDebrisTimeSteps = {0}".format(numDebrisTimeSteps)

    # Want to make the movie in all_points_delta_t.  How to get curSkyGrid into all_points_delta_t?
    # I think it used to be, but then I changed the code to make things more consistent.
    # Will have to write function that collapses probabilities into this timestep.
    # OR, simply pick one of the sub-timesteps to visualize
    # Prob at r_i in next x timesteps = 1 - (1 - probt1)(1 - probt2)...(1-probtx) just like consq probabilities?
    # OR just re-bin the curSkyGrid with the desired timestep.  Will it normalize out properly?
    if plotRisk:
        print "ERROR: Need to do some work make movies now that code has changed.  See comments above this error"
        sys.exit()


    # maxAcTimeStep = np.ceil((maxAircraftTrackTime - secondsFromMidnightUTC) / curMission['all_points_delta_t'])
    maxAcTimeStep = np.ceil((10.*60. / curMission['all_points_delta_t'])) # By 45 minutes, all AC have left the frame, use 50 tho

    figCounter = 1
    images = []
    for tx in range(numDebrisTimeSteps):

        # if not plotRisk:
        if tx >= maxAcTimeStep:
            # If you're not plotting the risk, then cut the animations once the AC data runs out.
            break

        if plotRisk:
            # Get the individual probability of impact grid
            probGrid = curSkyGrid.SendGridToPython(tx)
            if len(probGrid) == 0:
                print "No debris in NAS at time {0}".format(tx)
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

        if plotRisk:
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

        if plotRisk:
            # plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
            # cbarTicks = np.linspace(0, 1, 11, endpoint=True)
            # cbarTicks = np.log10(np.logspace(1e-5, 1, endpoint=True))
            cbarTicks = np.logspace(-4,0)

            # ax1.contourf(lonVec,latVec,probInterp, cbarTicks, cmap=plt.cm.jet, norm=LogNorm(), vmin=0., vmax=1.)

            # ax1.colorbar(ticks=cbarTicks)
            # plt.colorbar(ticks=cbarTicks, spacing='proportional')

        # print 'cbarTicks = ' + str(cbarTicks)
        plt.title("tx = {0:04d}, dt = {1} sec, altRange = [{2}, {3}]km".format(tx, curMission['all_points_delta_t'], lowAlt, highAlt))

        # Try to fix the axis
        # plt.axis([0, 6, 0, 20])
        plt.axis(corners)


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

        plt.savefig(curMission['GeneratedFilesFolder'] + "_tmp{0:04d}.png".format(figCounter))
        # images.append( (plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k'),) )
        figCounter += 1

        # plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
        # plt.contourf(lonVec,latVec,probInterp,15,cmap=plt.cm.jet)

    # im_ani = animation.ArtistAnimation(fig1, images, interval=1, repeat_delay=1, blit=True)

    # # im_ani = animation.Animation(fig1, images, blit=True)
    # # im_ani.save('GeneratedFiles/basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    # im_ani.save('doesntmatter.mp4', clear_temp=False)

    from subprocess import check_call
    print "Making Movie"
    check_call(["ffmpeg -y -r 5 -i {0}_tmp%04d.png -b:v 1800k {0}Contours.mp4".format(curMission['GeneratedFilesFolder'])], shell=True)
    # call(["ffmpeg -y -r 5 -i _tmp%04d.png -b:v 1800k ../../GeneratedFiles/Contours.mp4"], shell=True)
    print "Moving temporary images"
    check_call(["mv {0}_tmp*.png {0}MoviePngs/".format(curMission['GeneratedFilesFolder'])], shell=True)



    # plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
    # plt.savefig('GeneratedFiles/PlotGridProb')

    print 'exiting'
    sys.exit()















#     # ============ End of Risk Calculation =====================

#     plotRisk = True

#     '''
#     This is for making animations.  It's long and convoluted, sorry.
#     If you don't want one, you can stop the program here!
#     '''
#     if not makeAnimation:
#         sys.exit()

#     from scipy.interpolate import griddata
#     # import matplotlib
#     # matplotlib.use('Agg')  # Allows plot generation on server without X-windows
#     # import matplotlib.pyplot as plt

#     # import matplotlib.animation as animation
#     # import types # used in duck punching

#     np.set_printoptions(linewidth=120)

#     from matplotlib.colors import LogNorm
#     import matplotlib.gridspec as gridspec

#     # for plotting purposes
#     fig1 = plt.figure()


#     latShift = 0.086     # Increasing this number moves the test point downward
#     lonShift = -0.01      # Increasing this number moves the test point leftward
# #
# # # for plotting purposes
# # im = plt.imread('debrisField.png')
# # corners = [-98.2617 + lonShift, -92.1090 + lonShift, 29.6854 + latShift -0.04, 34.0093 + latShift]

#     debrisField = plt.imread('debrisField.png')
#     leftLon     = -98.2617 + lonShift
#     rightLon    = -92.1090 + lonShift
#     bottomLat   = 29.6854 + latShift -0.04
#     upperLat    = 34.0093 + latShift
#     corners = [leftLon,rightLon,bottomLat,upperLat]
#     # corners = [-98.2617,-92.1090,29.6854,34.0093]
#     # implot = ax1.imshow(debrisField,extent=corners)


#     # TIME LOOP WILL START HERE
#     altIX2plot  = 0
#     lowAlt      = altIX2plot * curMission['deltaZ']
#     highAlt     = (altIX2plot+1) * curMission['deltaZ']

#     if plotRisk:
#         # ASH the probabilities around
#         print "ASHING plotRisk"
#         curSkyGrid.generateASH(curMission['h1'], curMission['h2'])
#         print "DONE ASHING plotRisk"

#     # numDebrisTimeSteps will each be all_points_delta_t long
#     numDebrisTimeSteps = curSkyGrid.getNumRange()
#     print "numDebrisTimeSteps = {0}".format(numDebrisTimeSteps)

#     # maxAcTimeStep = np.ceil((maxAircraftTrackTime - secondsFromMidnightUTC) / curMission['all_points_delta_t'])
#     maxAcTimeStep = np.ceil((10.*60. / curMission['all_points_delta_t'])) # By 45 minutes, all AC have left the frame, use 50 tho

#     figCounter = 1
#     images = []
#     for tx in range(numDebrisTimeSteps):

#         # if not plotRisk:
#         if tx >= maxAcTimeStep:
#             # If you're not plotting the risk, then cut the animations once the AC data runs out.
#             break

#         if plotRisk:
#             # Get the individual probability of impact grid
#             probGrid = curSkyGrid.SendGridToPython(tx)
#             if len(probGrid) == 0:
#                 print "No debris in NAS at time {0}".format(tx)
#                 continue    # Debris hasn't made it to the NAS yet or there are no planes in the area of the debris

#             # The indices that correspond to the lowest alt level
#             altitudeLevels = np.unique(probGrid[:,0])
#             if len(altitudeLevels) < (altIX2plot+1):
#                 continue    # no probabilities to plot here.  move on
#             altIX = (probGrid[:,0] == altitudeLevels[altIX2plot])
#             altGrid = probGrid[altIX,:]

#             # Find the limits of the current grid in order to interpolate over them
#             minLat = np.min(altGrid[:,2])
#             maxLat = np.max(altGrid[:,2])
#             minLon = np.min(altGrid[:,1])
#             maxLon = np.max(altGrid[:,1])

#             # This is the mesh to interpolate over based on above limits
#             lonVec = np.linspace(minLon,maxLon,100)
#             latVec = np.linspace(minLat,maxLat,100)

#             # probInterp = griddata((altGrid[:,1], altGrid[:,2]), altGrid[:,3], (lonVec[None,:], latVec[:,None]), method='cubic')
#             probInterp = griddata((altGrid[:,1], altGrid[:,2]), altGrid[:,3], (lonVec[None,:], latVec[:,None]), method='linear')

#         print 'tx = {0} with a timestep of {1}'.format(tx, curMission['all_points_delta_t'])


#         # plt.clf()
#         # plt.imshow(debrisField,extent=corners)   # Hopefully this will wipe the old
#         # single_im = plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
#         # plt.text(0,1,"tx = {0:04d}".format(tx))

#         # # duck punching so that contour works
#         # def setvisible(self,vis):
#         #    for c in self.collections: c.set_visible(vis)
#         # single_im.set_visible = types.MethodType(setvisible,single_im,None)
#         # single_im.axes = plt.gca()
#         # images.append([single_im])

#         plt.clf()
#         gs = gridspec.GridSpec(1, 2,width_ratios=[15,1])

#         if plotRisk:
#             lowerBoundExponent = -4
#             upperBoundExponent = 0

#             # do the 2nd subplot, the pseudo colorbar, first
#             ax2 = plt.subplot(gs[1])
#             # np.logspace gives you logarithmically spaced levels -
#             # this, however, is not what you want in your colorbar
#             #
#             # you want equally spaced labels for each exponential group:
#             #
#             levls = np.array([])
#             for lix in range(lowerBoundExponent, upperBoundExponent):
#                 levls = np.concatenate((levls, np.linspace(10**lix,10**(lix+1),10)))

#             # levls = np.linspace(1,10,10)
#             # levls = np.concatenate((levls[:-1],np.linspace(10,100,10)))
#             # levls = np.concatenate((levls[:-1],np.linspace(100,1000,10)))
#             # levls = np.concatenate((levls[:-1],np.linspace(1000,10000,10)))

#             #
#             # simple x,y setup for a contourf plot to serve as colorbar
#             #
#             XC = [np.zeros(len(levls)), np.ones(len(levls))]
#             YC = [levls, levls]
#             CM = ax2.contourf(XC,YC,YC, levels=levls, norm = LogNorm())
#             # log y-scale
#             ax2.set_yscale('log')
#             # y-labels on the right
#             ax2.yaxis.tick_right()
#             # no x-ticks
#             ax2.set_xticks([])



#         ax1 = plt.subplot(gs[0])
#         ax1.imshow(debrisField,extent=corners)   # Hopefully this will wipe the old

#         if plotRisk:
#             # plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
#             # cbarTicks = np.linspace(0, 1, 11, endpoint=True)
#             # cbarTicks = np.log10(np.logspace(1e-5, 1, endpoint=True))
#             cbarTicks = np.logspace(-4,0)

#             # ax1.contourf(lonVec,latVec,probInterp, cbarTicks, cmap=plt.cm.jet, norm=LogNorm(), vmin=0., vmax=1.)

#             # ax1.colorbar(ticks=cbarTicks)
#             # plt.colorbar(ticks=cbarTicks, spacing='proportional')

#         # print 'cbarTicks = ' + str(cbarTicks)
#         plt.title("tx = {0:04d}, dt = {1} sec, altRange = [{2}, {3}]km".format(tx, curMission['all_points_delta_t'], lowAlt, highAlt))

#         # Try to fix the axis
#         # plt.axis([0, 6, 0, 20])
#         plt.axis(corners)


#         # Find the AC tracks at this time
#         for key in aircraftRecord:
#             curAC = np.array(aircraftRecord[key][0])
#             # goodIX = np.where((curAC[:,0] - secondsFromMidnightUTC) == 1.*(tx))
#             # print 'goodIX = ' + str(goodIX)

#             acType = aircraftRecord[key][1]

#             # goodIX = (curAC[:,0] - secondsFromMidnightUTC) == 1.*(tx)*mission1['all_points_delta_t']
#             timesHere       = range(tx*int(curMission['all_points_delta_t']), (tx+1)*int(curMission['all_points_delta_t']) )
#             goodIX          = np.in1d(curAC[:,0].ravel() - secondsFromMidnightUTC, timesHere).reshape(curAC[:,0].shape)
#             totalACLat      = curAC[goodIX,1]
#             totalACLon      = curAC[goodIX,2]
#             totalACLevel    = curAC[goodIX,3]

#             # print 'times here = ' + str(timesHere)
#             # print 'len = ' + str(len(totalACLat))
#             # Is there anybody here?
#             # if len(totalACLat) > 0:
#             for pt in range(len(totalACLat)):
#                 curLat = totalACLat[pt]
#                 curLon = totalACLon[pt]
#                 curAlt = totalACLevel[pt]
#                 # print '{0}, {1}'.format(curLon, curLat)
#                 if ((upperLat > curLat) & (curLat > bottomLat) & (leftLon < curLon) & (curLon < rightLon) & (lowAlt <= curAlt) & ( curAlt < highAlt)):
#                 #     print 'HERE!!!!'
#                     ax1.scatter(curLon ,curLat,marker='o',s=5, color='red')

#                 # If last pt
#                 if pt == (len(totalACLat)-1):
#                     ax1.annotate(acType, xy=(totalACLon[pt],totalACLat[pt]))

#         plt.savefig(curMission['GeneratedFilesFolder'] + "_tmp{0:04d}.png".format(figCounter))
#         # images.append( (plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k'),) )
#         figCounter += 1

#         # plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#         # plt.contourf(lonVec,latVec,probInterp,15,cmap=plt.cm.jet)

#     # im_ani = animation.ArtistAnimation(fig1, images, interval=1, repeat_delay=1, blit=True)

#     # # im_ani = animation.Animation(fig1, images, blit=True)
#     # # im_ani.save('GeneratedFiles/basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

#     # im_ani.save('doesntmatter.mp4', clear_temp=False)

#     from subprocess import check_call
#     print "Making Movie"
#     check_call(["ffmpeg -y -r 5 -i {0}_tmp%04d.png -b:v 1800k {0}Contours.mp4".format(curMission['GeneratedFilesFolder'])], shell=True)
#     # call(["ffmpeg -y -r 5 -i _tmp%04d.png -b:v 1800k ../../GeneratedFiles/Contours.mp4"], shell=True)
#     print "Moving temporary images"
#     check_call(["mv {0}_tmp*.png {0}MoviePngs/".format(curMission['GeneratedFilesFolder'])], shell=True)



#     # plt.contour(lonVec,latVec,probInterp,15,linewidths=0.5,colors='k')
#     # plt.savefig('GeneratedFiles/PlotGridProb')

#     print 'exiting'
#     sys.exit()




### Going to make some changes to SkyGrid, but these numbers shouldn't change!
# curMission['deltaXY']                   = .5    #km
# curMission['deltaZ']                    = NASkm/4.   #km
# curMission['h1']                        = 3.    # Smoothing parameters for the ASH.  Should be >= deltaXY
# curMission['h2']                        = 3.
# E120 = 0.000117191951219
# MD82 = 0.00161827291763
# B733 = 0.000948470675794
# B735 = 0.00107918905044
# MD90 = 0.0054834297398
# B738 = 0.0
# CRJ2 = 0.00273988944608
# CRJ2 = 0.00340425330452
# B733 = 0.000463863781848



















