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
vehicleName     = LaunchProviders.SS2
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

curMission['numNodes']                  = 1 # Will need to install pp to use more nodes
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
# curMission['launchAzimuth'] = 180.    #degrees, this is what it looks like on Google Earth
curMission['launchAzimuth'] = 0.    #degrees, Caveat: we're starting the launch at the point where SS2 is dropped
                                    #   so it's actually traveling BACK to the pad (180 degrees) but if I make this
                                    #   0 degrees, then I can rotate the envelope about the launch pad and get the angles
                                    #   to look right.  It's kind of a hack, but only kinda.


#
# curFileName = '/Users/marian/Documents/Research/Prop3Dof/CythonFiles/OtherPythonFiles/FootprintLibrary/workInProgressTJC01/SS2_America_workInProgressHEALTH.dat'
# curFootprint = ceb.PyFootprint(curFileName, True)
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
    t_hi = 120. #420.

    TJC.MonteCarlo_until_tfail(curMission, profiles, t_lo, t_hi)

# # ## Find the time until the airspace can become reactive
# minTime = 60.
# maxTime = 180.
# tProactive = TJC.FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime)
# print "tProactive = {0}\n".format(tProactive)
# TJC.PlotNominalTrajectories(profiles, curMission, maxTime)
# sys.exit()

if not debug:
    footprintStart = 0.
    # footprintUntil = 200.
    footprintUntil = 420.
    footprintIntervals = 60.

    footprintTotal = TJC.GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals)
    vehicleNotes = vehicleNotes + 'HealthFlash' + str(int(footprintIntervals))

    # footprintTotal = TJC.GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals)
    # vehicleNotes = vehicleNotes + 'NoHealth' + str(int(footprintIntervals))

    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes)

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(outfileStr)
    sys.exit()


else:
    print 'DEBUGGING NEW WAY OF GENERATING ENVELOPES'
    sys.exit()

    # EVstrike, Fprint = TJC.makeFootprintFromTimes(curMission, tfailSec, tfailSec)
    # Let's start by looking at the first 10 seconds for the very first fail time
    # Load the precomputed footprint


    ''' Prototype Flash'''
    footprintStart = 0.
    # footprintUntil = 50.
    footprintUntil = 420.
    footprintIntervals = 60.
    footprintTotal = []

    for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
        timelo = footprintStart + ix*footprintIntervals
        timehi = np.min( (footprintStart + (ix+1)*footprintIntervals, footprintUntil) )

        print 'TIMES: From {0} to {1}'.format(timelo, timehi)
        EVstrike, curFootPrint = TJC.makeFootprintFromTimes(curMission, timelo, timehi)
        print 'EV =  ' + str(EVstrike)

        # Now take that footprint and...
        # Smooth it out to a single timestep
        numRange = curFootPrint.getNumRange()
        curFootPrint.SmoothedOut(numRange)  # This will make footprintDelaT = numRange, and then change numRange to = 1

        numRange = curFootPrint.getNumRange()
        FPDeltaT = curFootPrint.getDeltaT()

        # Resize the deltaT to be only the length of the interval
        #   So if we're making an envelope at each second, then the footprint should be chopped at 1 second
        #   If we're combining times, like every 5 seconds or every minute, then it should be 5s or 60s
        curFootPrint.ChopTimeAt(footprintIntervals)

        # Translate the footprint forward to tfailSec
        if timelo > 0:
            curFootPrint.SlideFootprintBySeconds(timelo)

        # Merge it with the others
        if ix == 0:
            footprintTotal = curFootPrint
        else:
            print '\n\nMERGE'
            footprintTotal.MergeFootprintVectors(curFootPrint)


        # # Print to GE
        # debugFolder = 'GeneratedFiles/Sandbox/'
        # vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, timelo, timehi)
        # curFootPrint.ExportGoogleEarth(debugFolder + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)

        # # Fprint.SmoothedOut(footprintIntervals)
        #
        # # Fprint.SmoothedOut()
        # curFootPrint.ExportGoogleEarth('GeneratedFiles/PythonGE_' + str(timelo) + 'To'
        #                                       + str(timehi) + 'FootprintSMOOTH.kml', yyyy, mm, dd, hour, min)

    # Just to be safe(?), set the params we need in order to translate / rotate
    footprintTotal.SetAzimuthDeg(curMission['launchAzimuth'])
    footprintTotal.SetLaunchLatDeg(curMission['launchLat'])
    footprintTotal.SetLaunchLonDeg(curMission['launchLon'])

    # Print to GE
    vehicleFileName = '{0}_{1}_{2}'.format(vehicleName, launchLocation, vehicleNotes + '_HEALTH')

    footprintTotal.ExportGoogleEarth(curMission['footprintLibrary'] + vehicleFileName + '.kml', yyyy, mm, dd, hour, min)
    outfileStr = curMission['footprintLibrary'] + vehicleFileName + '.dat'
    footprintTotal.StoreFootprintAsVector(outfileStr)



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



































'''
#
This file creates an envelope for a Spaceship2 launch and reentry
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



# ============= Define the Mission ====================== #
# from ReadMission import readInput

# Initialize the mission
propagationParamFile = []
precomputedParamFile = 'HTHLwCarrier.txt'
pathToMissionFiles = 'Files/SS2/'

# Planet info
omegaE = 7.2921158494529352e-05     # rad/s
planetModel = 0

mission1 = dict(propagationParamFile = propagationParamFile, precomputedParamFile = precomputedParamFile,
                pathToMissionFiles = pathToMissionFiles, omegaE = omegaE, planetModel = planetModel)
mission1 = TJC.InitializeMission(mission1)

# Decide if you want to implement this the same way as for vehicles we have thrust profiles for
#mission1 = readInput('nominalParam_new.txt', pathToMissionFiles)
# mission1['pathToMissionFiles'] = pathToMissionFiles

# Generate a realistic profile
#from Files.SS2.failProfile import failProfile, failProfileSeconds   # This should go in the readInput file
#mission1['failProfile'] = failProfile
#mission1['failProfileSeconds'] = failProfileSeconds


mission1['pFail'] = 0.01     # Probability that vehicle will fail somewhere
mission1['reactionTimeMinutes'] = 5

# My Stuff
#mission1 = dict()
mission1['GeneratedFilesFolder'] = pathToMissionFiles + 'GeneratedFiles/'
GeneratedFilesFolder = mission1['GeneratedFilesFolder']

#mission1['LrhcFile'] = GeneratedFilesFolder + 'CapeLRHC.dat'
mission1['GoogleEarthFile'] = GeneratedFilesFolder + 'PythonGE.kml'

mission1['deltaT'] = 1.    #Seconds, this is the time resolution of a propagated trajectory
mission1['deltaTFail'] = 5.0   #Seconds, this is how often we explode the rocket

mission1['debrisPickleFolder'] = GeneratedFilesFolder + 'debrisPickleFolder'
mission1['footprintVectorFolder'] = GeneratedFilesFolder + 'footprintVectorFolder'

# Make sure that the directory for holding the general Generated files exists
folderPath = os.path.abspath(GeneratedFilesFolder)
if not os.path.exists(folderPath):
    os.makedirs(folderPath)

mission1['exportGE'] = True
mission1['exportFACET'] = False

# Frisco Stuff
mission1['debrisCatPath'] = pathToMissionFiles + 'DebrisCatalog/'
mission1['debrisCatFile'] = 'SI_Halcon9_2nd.dat'
mission1['atmosphere'] = friscoFiles + 'AtmoProfiles/special.txt'

# Updated, but not sure i need these
mission1['launchLat'] = 32.989948
mission1['launchLon'] = -106.975451
mission1['launchAlt'] = 15.24   #km

# Not updated, also not sure i need these
mission1['initialUTC'] = 156.84861111111113
mission1['launchAzimuth'] = 0.976192685683


# Changing this to seconds!!!  Make sure that this timestep GREATER THAN OR EQUAL to deltaT, right?
mission1['all_points_delta_t'] = 1.0    #seconds, this will be the time resolution of a compact envelope



# KDE / Footprint Parameters
mission1['deltaXY'] = 5.   #km
mission1['deltaZ'] = 20.   #km

mission1['thresh'] = 1e-7
mission1['h1'] = 30.
mission1['h2'] = 30.


yyyy = 2014
mm = 1
dd = 1
hour = 2
min = 13

mission1['ExportDate'] = [yyyy, mm, dd, hour, min]

debrisPickleFolder = mission1['debrisPickleFolder']
deltaTFail = mission1['deltaTFail']

mission1['numPiecesPerSample'] = 5

# ======================= Begin Computations ============================ #



########### Fold that into existing infrastructure

# Precompute some Atmosphere and Trajectory profiles
fresh = True

profiles = []
if (fresh):
    # Should really move all the important mission stuff into this if-statement and wrap it up into the montecarlo dictionary
    
    numTrajSamples = 1
    numWindSamples = 5
    
    # I only need to generate wind profiles here, since i'm not going to worry about multiple nominal trajectories yet
    # Could / should probably anticipate doing it though andjust replicate the single trajectory here to conform with the existing infrastrcture
    
    atmStorage, stateVecStorage, thetagStorage, tfailStorage = TJC.GenerateWindTrajProfiles(mission1, numTrajSamples, numWindSamples)
    profiles = dict(atmStorage = atmStorage, stateVecStorage = stateVecStorage, thetagStorage = thetagStorage, tfailStorage = tfailStorage,
                    numTrajSamples = numTrajSamples, numWindSamples = numWindSamples)


    # Test Velocity Profile
    whichTraj = 0
    whichWind = 0
    X = profiles['stateVecStorage'][whichWind][whichTraj][0]
    Y = profiles['stateVecStorage'][whichWind][whichTraj][1]
    Z = profiles['stateVecStorage'][whichWind][whichTraj][2]

    Vx = profiles['stateVecStorage'][whichWind][whichTraj][3]
    Vy = profiles['stateVecStorage'][whichWind][whichTraj][4]
    Vz = profiles['stateVecStorage'][whichWind][whichTraj][5]

    VrelMag = profiles['stateVecStorage'][whichWind][whichTraj][6]

    Rnorm = np.sqrt(X*X + Y*Y + Z*Z)*(1e-3) - (6371)  #rEarth in km
    Vnorm = np.sqrt(Vx*Vx + Vy*Vy + Vz*Vz)*(1e-3)

    print np.max(Rnorm)
    print np.max(Vnorm)

    # Veci = np.array([378.34173789620399, -14.650954707525671, -101.34663177561015])
    # Vrotfinal = np.array([388.10047169840857, 34.57221526682423, 0.0])

    # print np.linalg.norm(Veci - Vrotfinal)

    sys.exit()

    import pickle
    output = open('localProfiles.pkl', 'wb')
    pickle.dump(profiles,output)
    output.close()
    
    tfail = 135.0
    
    TJC.MonteCarlo_until_tfail(mission1, profiles, tfail)
    

    

else:
    import pickle
    profiles = pickle.load(open('localProfiles.pkl','rb'))







