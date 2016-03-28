import os
import sys

# Points to the python scripts needed from Francisco
# friscoFiles = '../../../Prop3Dof/FriscoDebris/pythonFiles/'

# curPath = os.path.dirname(__file__) + "/"
# friscoFiles = 'FriscoLegacy/'
# sys.path.append(curPath + friscoFiles)

# friscoDebris = os.path.abspath(curPath + "../../build/")
# sys.path.append(friscoDebris)


# print sys.path



import numpy as np
import CompactEnvelopeBuilder as ceb
# import getPropTraj as traj
from FriscoLegacy import AtmosProfile as AP
from FriscoLegacy import debrisReader as DR
from FriscoLegacy import debrisPropagation as dp



# from scipy.io import loadmat
from FriscoLegacy import orbitTools
import pickle

# import pdb

import math

# TODO: Make sure the cell volume is SMALLER than Wilde's assumed aircraft density, otherwise P_aircraft > 1

# TODO: IMPORTANT: Make it such that the top of the zBin is ALWAYS the NAS ceiling.  If we allow the highest bins to technically
#   stretch above the NAS, that makes their volume larger than the volume over which probabilities were calculated.
#   Make all bins equal height from zero to NAS.  Similarly, stop using numZBins in C++ and just use the height so that
#   you're always working with true volumes.  That seems much less confusing.

# TODO: Should have option to sample new wind profile for every debris piece to improve convergence
# I think i did it this way before so that debugging would be easier, but now i think it's suboptimal
#   because it keeps us from properly capturing the variation in the wind in the monte carlo sim

# TODO: Make sure that the three important time resolutions nest well
# TODO: Make sure that the grid spacing and the smoothing parameters nest
# TODO: If I change all_points_delta_t to be other values, i get a segfault
# TODO: I commented out the line in skygrid that forces deltaZ <= NASkm.  Handle this issue properly

# TODO: If using a reactionTime, then no need to propagate debris beyond that time.  Could save TONS of time and keep
#   me from hitting those weird memory errors on Zion.




'''
Defines the file structure for folders / files that are specific to this vehicle.
In theory, all vehicles will have the same file structure and this section would
    be the same across all main files.  Also, should bring pathToMissionFiles
    into this section
'''
def SetupOutputFolders(curMission, tempDir, outputDir, vehicleName, launchLocation):

    # These hold output files, create them in a moment if they don't already exist
    curMission['GeneratedFilesFolder']    = tempDir + "{0}_{1}/".format(vehicleName, launchLocation)
    curMission['debrisPickleFolder']      = curMission['GeneratedFilesFolder']  + 'debrisPickleFolder'
    curMission['footprintVectorFolder']   = curMission['GeneratedFilesFolder']  + 'footprintVectorFolder'
    curMission['facetFolder']             = curMission['GeneratedFilesFolder']  + 'facetFolder/'
    curMission['footprintLibrary']        = outputDir  + 'footprintLibrary/'

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

    # del folderPath      # Just to be safe!

# Values taken from frisco's orbitTools.py
def geodetic2spherical(gdlat,lon,hs):
    ECEF = orbitTools.latlonalt2ECEF(gdlat,lon,hs,1)
    return orbitTools.ECEF2latlonalt(ECEF,0)


def createAtmoPickle(atmoFolder, atmoFile, pickleName):
# atmoFile = 'special.txt'

    # getting atmospheric profile
    altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,\
        densitySDList,uSDList,vSDList,wSDList,nPoints = AP.readGramAtmos(atmoFolder + atmoFile)

    # making sure atmospheric profile is decreasing in altitude
    altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,\
        densitySDList,uSDList,vSDList,wSDList,nlist \
        = AP.maxAltitudeAtmosFix(altitudeList,densityMeanList,uMeanList,vMeanList,\
                                 wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints)

    # [altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist] \
    #     = pickle.load(open(curMission['atmospherePickle'],'rb'))

    output = open(atmoFolder + pickleName, 'wb')
    pickle.dump([altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist],output,2)
    output.close()




# def testASH(mission1, deltaTFail, timelo, timehi):
#     # Let's make the updating footprint for oneToTwo minutes
#     #    timehi = 115.
#     #    timelo = 60.
#
#     debrisPickleFolder = mission1['debrisPickleFolder']
#
#     numGridsHere = int(np.round((timehi - timelo)/deltaTFail) + 1)
#
#     for ix in range(numGridsHere):
#         tfailSec = timelo + ix*deltaTFail
#         print 'tfailSec = ' + str(tfailSec)
#
#         input = open(debrisPickleFolder + '/mpc_' + str(tfailSec) + '.pkl', 'rb')
#         cur_mpc = pickle.load(input)
#         input.close()
#
#         curPointCloud = ceb.PyPointCloud(cur_mpc, tfailSec)
#         curSkyGrid    = ceb.PySkyGrid(curPointCloud, mission1['deltaXY'], mission1['deltaXY'], 5)
#
#
#         curSkyGrid.ASH2()
#
# #        if (ix == 0):
# #            oneToTwoFootprint = ceb.PyFootprint(curSkyGrid)     # create the footprint
# #        else:
# #            myFootprint = ceb.PyFootprint(curSkyGrid)
# #            oneToTwoFootprint.MergeFootprintVectors(myFootprint)        # merge the new footprint into the old one
#
# #    oneToTwoFootprint.ProjectAllPointsDown()
# #    return oneToTwoFootprint
















'''
GenerateEnvelopes_HealthFlash: Calculates the envelopes to a threshold with a reaction time and then chains them all
    together.  All individual envelopes are projected down into a single timestep which is the length of
    footprintIntervals.  These footprints, when chained together, will simulate instantaneous health monitoring.
    Example: If the vehicle passes through t=2 nominally, then E_{t=2} will be turned off at time t=3.
'''
def GenerateEnvelopes_HealthFlash(curMission, footprintStart, footprintUntil, footprintIntervals):

    footprintTotal = []

    for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
        timelo = footprintStart + ix*footprintIntervals
        timehi = np.min( (footprintStart + (ix+1)*footprintIntervals, footprintUntil) )

        print 'TIMES: From {0} to {1}'.format(timelo, timehi)
        EVstrike, curFootPrint = makeFootprintFromTimes_InstantaneousOnly(curMission, timelo, timehi)
        print 'EV =  ' + str(EVstrike)

        # Now take that footprint and...
        # Smooth it out to a single timestep

        # This is what the footprint will look like if there was an "error" in generating it
        if curFootPrint == []:
            print "\n\n===========\nEMPTY FOOTPRINT!!!\n==============\n\n"
            continue

        numRange = curFootPrint.getNumRange()
        curFootPrint.SmoothedOut(numRange * curMission['all_points_delta_t'])  # This will make footprintDelaT = numRange, and then change numRange to = 1

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

    footprintTotal.SmoothedOut(0)   #I believe this will simply smooth the footprints and not alter the timesteps

    # Just to be safe(?), set the params we need in order to translate / rotate
    footprintTotal.SetAzimuthDeg(curMission['launchAzimuth'])
    footprintTotal.SetLaunchLatDeg(curMission['launchLat'])
    footprintTotal.SetLaunchLonDeg(curMission['launchLon'])

    return footprintTotal





'''
GenerateEnvelopes_NoHealth: Calculates the envelopes to a threshold with a reaction time and then chains them all
    together.  If the vehicle passed through a timestep without exploding, this architecture doesn't care.
'''
def GenerateEnvelopes_NoHealth(curMission, footprintStart, footprintUntil, footprintIntervals):
    footprintTotal = []

    for ix in range(int(np.ceil((footprintUntil-footprintStart)/footprintIntervals))):
        timelo = footprintStart + ix*footprintIntervals
        timehi = np.min( (footprintStart + (ix+1)*footprintIntervals - curMission['deltaTFail'], footprintUntil) )
        EVstrike, Fprint = TJC.makeFootprintFromTimes_InstantaneousOnly(curMission, timelo, timehi)

        print 'EV =  ' + str(EVstrike)

        Fprint.SmoothedOut(footprintIntervals)
        # Fprint.SmoothedOut()
        # Fprint.ExportGoogleEarth('GeneratedFiles/PythonGE_' + str(timelo) + 'To'
        #                                       + str(timehi) + 'FootprintSMOOTH.kml', yyyy, mm, dd, hour, min)

        if ix == 0:
            footprintTotal = Fprint
        else:
            footprintTotal.MergeFootprintVectors(Fprint)


    # Just to be safe(?), set the params we need in order to translate / rotate
    footprintTotal.SetAzimuthDeg(curMission['launchAzimuth'])
    footprintTotal.SetLaunchLatDeg(curMission['launchLat'])
    footprintTotal.SetLaunchLonDeg(curMission['launchLon'])

    return footprintTotal




def makeFootprintFromTimes_InstantaneousOnly(mission1, timelo, timehi):
    # Let's make a footprint that bounds all debris generated between lo and hi with a boundary that's <= thresh

    # Unpack some things
    GeneratedFilesFolder        = mission1['GeneratedFilesFolder']
    failProfile                 = mission1['failProfile']
    failProfileSeconds          = mission1['failProfileSeconds']
    [yyyy, mm, dd, hour, min]   = mission1['ExportDate']
    deltaTFail                  = mission1['deltaTFail']
    pFail = mission1['pFail']

    # Will take a min on this later, so it will get set to the max nodes available
    numNodesEnvelopes = 100
    if mission1.has_key('numNodesEnvelopes'):
        # Unless we're capping the nodes at something even less
        numNodesEnvelopes = mission1['numNodesEnvelopes']

    # Incoming times should ideally nest with deltaTFail 
    # Will form the envelope that corresponds to failures [timelo, timehi)
    # e.g. [100, 105) = failures from 100,101,102,103,104
    if (timehi - timelo) % deltaTFail != 0:
        print "ERROR: Make your incoming timesteps nest!"
        sys.exit()

    # Find out how many debris data files there are for these times
    # Already enforced their nesting with mod, so this should be safe
    numGridsHere = int(np.round((timehi - timelo)/deltaTFail)) # TODO: Fix this when removing overlapping times


    # # Actually, this method isn't valid for anything other than instantaneous, so make sure!
    # if numGridsHere > 1:
    #     print "ERROR: This function is only valid for instantaneous!  No merging timesteps allowed"
    #     sys.exit()
    # It is okay to merge timesteps, say if we're getting health updates every five seconds but the compact
    #  envelope *nominally* evolves every 60 seconds.

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

    # genFootprint should not be worried about pFail!!!

    jobs = []   # Holds the results of the envelope creation
    if numNodesEnvelopes == 1:
        # This means we're running into memory or malloc issues with pp, so don't use it
        jobs = [genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) for ix in range(len(timeRange))]
    else:
        # Use the specified number of nodes in parallel
        # Now that we know which times to propagate, do them in parallel
        import pp
        job_server = pp.Server()
        numPossibleNodes = job_server.get_ncpus()
        job_server.set_ncpus(np.min((numPossibleNodes, numNodesEnvelopes)))
        print 'number of cpus = ' + str(job_server.get_ncpus())

        # Create all of the footprints in parallel
        print 'Youre about to submit jobs in parallel.  BE AWARE that cout statements will cause python to malloc error here.  Dont know why.'
        print 'So if you get a malloc error, be suspicious that you hit an error which triggered a cout'
        jobs = [job_server.submit(genFootprint, \
                                  args=(mission1, timeRange[ix], pFailThisTimestepVec[ix]), \
                                  callback=finished002 \
                                  ) for ix in range(len(timeRange))]
        job_server.wait()
        job_server.print_stats()


    vals        = []    # Holds the EV_Strike values
    fileNames   = []    # Holds the filenames of the footprint vector .dat files
    if numNodesEnvelopes == 1:
        # This is a list of tuples
        for job in jobs:
            print job
            vals.append(job[0])
            fileNames.append(job[1])
    else:
        # This is a function returned from pp
        for job in jobs:
            print job()
            vals.append(job()[0])
            fileNames.append(job()[1])
        #job_server.destroy()    #Hopefully, this will close the server and we won't get errors about too many open files

    try:
        # Need to open up all the footprints and merge them
        print 'Merging Footprints'
        totalFootPrint = ceb.PyFootprint(fileNames[0], True)
        print fileNames[0]

        for curFile in fileNames[1:]:
            curFootPrint = ceb.PyFootprint(curFile, True)
            totalFootPrint.MergeFootprintVectors(curFootPrint)  # merge the new footprint into the old one
                                                                #   The points will need to be smoothed eventually

    except IndexError:
        if (len(timeRange) == 0):
            # There were no non-zero probability events
            vals = 0.
            totalFootPrint = []
        else:
            raise
    except:
        raise

    # totalFootPrint.ProjectAllPointsDown()

    return np.max(vals), totalFootPrint






# def makeFootprintFromTimes(mission1, timelo, timehi):
#     # Let's make a footprint that bounds all debris generated between lo and hi with a boundary that's <= thresh

#     # Explosions at tfailSec = 0 is NOT allowed.  Must be strictly > zero.
#     # For now, I'm saying that with deltaTFail = 5, a tfailSec = 20 represents the debris
#     #   generated during the statetime period from t = 15 UP TO 20sec.
    
#     # Unpack some things
#     # debrisPickleFolder          = mission1['debrisPickleFolder']
#     # footprintVectorFolder       = mission1['footprintVectorFolder']
#     GeneratedFilesFolder        = mission1['GeneratedFilesFolder']
#     failProfile                 = mission1['failProfile']
#     failProfileSeconds          = mission1['failProfileSeconds']
#     # deltaXY                     = mission1['deltaXY']
#     # deltaZ                      = mission1['deltaZ']
#     [yyyy, mm, dd, hour, min]   = mission1['ExportDate']
#     # ExportDate                  = mission1['ExportDateDT']

#     deltaTFail                  = mission1['deltaTFail']
#     # thresh                      = mission1['thresh']
#     # h1                          = mission1['h1']
#     # h2                          = mission1['h2']

#     # Will take a min on this later, so it will get set to the max nodes available
#     numNodesEnvelopes = 100
#     if mission1.has_key('numNodesEnvelopes'):
#         # Unless we're capping the nodes at something even less
#         numNodesEnvelopes = mission1['numNodesEnvelopes']
    
#     pFail = mission1['pFail']
#     # reactionTimeSeconds = mission1['reactionTimeSeconds']
    
#     numGridsHere = int(np.round((timehi - timelo)/deltaTFail) + 1)

#     # Figure out which failure times are worth propagating (i.e. they have a nonzero probability of happening)
#     timeRange = []
#     pFailThisTimestepVec = []
#     lastIndices = []
#     for ix in range(numGridsHere):
#         tfailSec = timelo + ix*deltaTFail
        
#         # failProfileSeconds likely has a smaller timestep than deltaTFail, so will need find the failProfileSeconds
#         #  that fall within the range [tfailSec - deltaTFail, tfailSec).  This is using the probability of all the times
#         #  previous to tfailSec (by deltaTFail) up to but-not-including tfailSec to approximate the probability of failing at tfailSec
#         indices = np.where((failProfileSeconds >= (tfailSec - deltaTFail)) & (failProfileSeconds < tfailSec))[0]
#         pFailThisTimestep = np.sum(failProfile[indices])

#         # If pFailThisTimestep < thresh, then the expected number of strikes for this event will surely be less than thresh
#         if pFailThisTimestep > 0.:
#             timeRange.append(tfailSec)
#             pFailThisTimestepVec.append(pFailThisTimestep * pFail)
#             lastIndices = indices

#     print 'times are ' + str(timeRange)
#     print 'probs are ' + str(pFailThisTimestepVec)
#     print 'lastIndices are ' + str(lastIndices)


#     # Make sure you can do one
#     # TODO: Don't even bother with this if you're only using one node
#     print 'DOING A TEST RUN BEFORE PARALLEL'
#     try:
#         testIX = -1     # Look at the highest time

#         print 'timeRange = {0}'.format(timeRange)
#         print '{0} - {1}, {2}'.format(timelo, timehi, numGridsHere)

#         # Make sure that there are no cout statements that'll trip up the parallel methods
#         # print '     THERE SHOULD BE NOTHING BETWEEN THIS LINE'
#         # EV_strike, outfileStr = genFootprint(mission1, timeRange[testIX], pFailThisTimestepVec[testIX])
#         # print '     AND THIS LINE'



#         # # Print out the footprint generated by this test run
#         # GEFileName = GeneratedFilesFolder + 'PythonGE_' + str(timeRange[testIX]) + '.kml'
#         # debugFootprint = ceb.PyFootprint(outfileStr, True)
#         # debugFootprint.ExportGoogleEarth(GEFileName, yyyy, mm, dd, hour, min)
#         # print 'GEfile = ' + GEFileName
#     except:
#         if len(timeRange) == 0:
#             # The probabilty of anything happening here is zero
#             EV_strike = 0.
#         else:
#             raise
#     # print 'EV_strike = ' + str(EV_strike)


#     jobs = []   # Holds the results of the envelope creation
#     if numNodesEnvelopes == 1:
#         # This means we're running into memory issues with pp, so don't use it
#         jobs = [genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) for ix in range(len(timeRange))]
#     else:
#         # Use the specified number of nodes in parallel
#         # Now that we know which times to propagate, do them in parallel
#         import pp
#         job_server = pp.Server()
#         numPossibleNodes = job_server.get_ncpus()
#         job_server.set_ncpus(np.min((numPossibleNodes, numNodesEnvelopes)))
#         print 'number of cpus = ' + str(job_server.get_ncpus())

#         # Create all of the footprints in parallel
#         print 'Youre about to submit jobs in parallel.  BE AWARE that cout statements will cause python to malloc error here.  Dont know why.'
#         print 'So if you get a malloc error, be suspicious that you hit an error which triggered a cout'
#         jobs = [job_server.submit(genFootprint, \
#                                   args=(mission1, timeRange[ix], pFailThisTimestepVec[ix]), \
#                                   callback=finished002 \
#                                   ) for ix in range(len(timeRange))]
#         job_server.wait()
#         job_server.print_stats()
#         job_server.destroy()    #Hopefully, this will close the server and we won't get errors about too many open files


#     vals        = []    # Holds the EV_Strike values
#     fileNames   = []    # Holds the filenames of the footprint vector .dat files
#     if numNodesEnvelopes == 1:
#         # This is a list of tuples
#         for job in jobs:
#             print job
#             vals.append(job[0])
#             fileNames.append(job[1])
#     else:
#         # This is a function returned from pp
#         for job in jobs:
#             print job()
#             vals.append(job()[0])
#             fileNames.append(job()[1])


#     try:
#         # Need to open up all the footprints and merge them
#         print 'Merging Footprints'
#         totalFootPrint = ceb.PyFootprint(fileNames[0], True)
#         print fileNames[0]

#         for curFile in fileNames[1:]:
#             curFootPrint = ceb.PyFootprint(curFile, True)
#             totalFootPrint.MergeFootprintVectors(curFootPrint)  # merge the new footprint into the old one
#                                                                 #   The points will need to be smoothed eventually

#     except IndexError:
#         if (len(timeRange) == 0):
#             # There were no non-zero probability events
#             vals = 0.
#             totalFootPrint = []
#         else:
#             raise
#     except:
#         raise

#     # totalFootPrint.ProjectAllPointsDown()

#     return np.max(vals), totalFootPrint


def finished002(val):
#    print "done with something"
    print "done with " + val[1]


def genFootprint(curMission, tfailSec, curPFail):
    """For a non-reentry.  This takes a time of failure and its associated probability and generates a 
        a footprint around the debris that would be created from this failure.  Writes footprint to
        file and returns the danger metric and the name of the file."""

    print "tfailSec {0}, curPFail {1}".format(tfailSec, curPFail)

    ExportDateDT            = curMission['ExportDateDT']
    # reactionTimeSeconds     = curMission['reactionTimeSeconds']       # Where did this go?
    h1                      = curMission['h1']
    h2                      = curMission['h2']
    debrisPickleFolder      = curMission['debrisPickleFolder']
    footprintVectorFolder   = curMission['footprintVectorFolder']
    thresh                  = curMission['thresh']
    cumulative              = curMission['cumulative']
    whichProbability        = curMission['whichProbability']


    from CompactEnvelopeBuilder import PySkyGrid, PyPointCloud, PyFootprint
    curSkyGrid              = []     # Scope these here
    arefMeanList            = []
    numberOfPiecesMeanList  = []

    # Open up the debris
    if curMission.has_key('isReentry') and curMission['isReentry'] == True:
        prevAref = []
        prevMean = []
        counter = 0

        for windIX in range(curMission['numWindProfiles']):
            inFileName = '{0}/mpc_{1}_{2}.pkl'.format(curMission['debrisPickleFolder'], windIX, int(tfailSec))
            # print inFileName

            input = open(inFileName, 'rb')
            cur_mpc = pickle.load(input)
            input.close()

            arefMeanList = cur_mpc['arefMeanList']
            numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']

            # Package them up into a PointCLoud
            # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
            curPointCloud = PyPointCloud(cur_mpc, tfailSec, curMission)

            if counter == 0:
                prevAref = arefMeanList
                prevMean = numberOfPiecesMeanList

                # Place the cloud into a fresh Grid
                curSkyGrid    = PySkyGrid(curMission, curPointCloud)
            else:
                if (arefMeanList != prevAref):
                    print 'ERROR: Using different area lists'
                    # Comments in the montecarlo functions indicate that this is no longer needed and is []
                    sys.exit()
                elif (numberOfPiecesMeanList != prevMean).all():
                    print 'ERROR: Using different mean lists'
                    sys.exit()
                else:
                    #Everything is fine, append to existing grid
                    1
                    curSkyGrid.IncorporatePointCloudIntoGrid(curPointCloud)

            counter = counter + 1

    else:
        inFileName = '{0}/mpc_{1}.pkl'.format(debrisPickleFolder, str(tfailSec))
        input = open(inFileName, 'rb')
        cur_mpc = pickle.load(input)
        input.close()

        arefMeanList = cur_mpc['arefMeanList']
        numberOfPiecesMeanList = cur_mpc['numberOfPiecesMeanList']

        # This is [total number of pieces simulated within this mpc] / [number of debris categories in this mpc]
        # TODO: If all_points_delta_t != debrisDeltaT, then we'll be double-counting here.
        # numDebrisPerIXSimulated = cur_mpc['numPieces']/len(numberOfPiecesMeanList)

        # Package them up into a PointCLoud
        # NOTE!!!  Inside the PointCloud constructor we apply the reactionTime which is NO LONGER HARDCODED!!!
        curPointCloud = PyPointCloud(cur_mpc, tfailSec, curMission)

        # Place the cloud into a Grid
        curSkyGrid    = PySkyGrid(curMission, curPointCloud)

    # Do the ASH
    # Must do the whole thing up-front.  On the fly only works with risk calculations at certain predetermined points.
    print 'ASHING'
    curSkyGrid.generateASH(h1, h2)

    # print 'tfailsec = ' + str(tfailSec)
    # return -1, './GeneratedFiles/footprintVectorFolder/fpVec_20.0.dat'

    # TODO: Move this ac density map stuff to a graveyard.
    # # Get the lat/lons of the filled cells
    # print 'createEmptyAircraftDensityMap'
    # latlonArray = curSkyGrid.createEmptyAircraftDensityMap()

    # Using density maps are just not a good idea.  
    # if useAircraftDensityMap:
    #     # With those lat/lons, find the probAircraft for each cell
    #     from AircraftDensityMap import AircraftDensityMap as ADM
    #     density = ADM()
    #     densityArray = density.getDensity(latlonArray, ExportDateDT)
    #     print 'THIS IS PROBABLY A DENSITY AND NOT A PROBABILITY.  FIX THIS!!!'
    #     sys.exit()

    #     # Send that information back into C++
    #     curSkyGrid.populateAircraftDensityMap(densityArray, len(densityArray))
    # else:
    #     fourNM2 = 13.72                     #// 4 (n.m.)^2 * (1.852 km/nm)^2 = 13.72 km^2
    #     aircraftDensity = 1./fourNM2        #// [prob/km^2] Paul Wilde's assumed aircraft density (1 every 4nm^2)
    #     cellArea = deltaXY*deltaXY
    #     probOfAirplaneInCell = aircraftDensity * cellArea;

    #     import numpy as np
    #     # Send that information back into C++
    #     print 'populateAircraftDensityMap'
    #     curSkyGrid.populateAircraftDensityMap(np.array([probOfAirplaneInCell]), -1)

    # After uploading the density map, we have to generate the hazard probabilities
    print 'generateHazardProbabilities'
    curSkyGrid.generateHazardProbabilities(numberOfPiecesMeanList)


    # Now apply the current definition of 'cumulative' for the chosen type of probability
    if cumulative == 'FAA':
        # # This will require coarsening the grid
        # newDeltaXY   = 3.5      #//[km]
        # newDeltaZ    = 20.      #//[km]  This is higher than NASkm, but I need the values to nest, so hopefully this is fine
        #
        # print 'SHIT GUYS!
        # raise RuntimeError
        newDeltaXY   = -1      #//[km]
        newDeltaZ    = -1

        # Inside this function, the grid will be permanently coarsened.  I don't think it's a good idea to try
        #   to alter the curMission here because if it's shared memory than that could potentially mess up other
        #   calculations that are happening in parallel threads.  Wait until whole thing is over and then do it.
        print 'generateAllPoints_CumulativeFAA'
        EV_strike = curSkyGrid.generateAllPoints_CumulativeFAA(thresh, whichProbability, curPFail)

    elif cumulative == 'TJC':
        # EV_strike = curSkyGrid.generateAllPoints_CumulativeTJC(thresh, whichProbability)
        print "ERROR: TJC Cumulative no longer exists"
        sys.exit()



    # print 'EV_Strike = ' + str(EV_strike)
    outfileStr = footprintVectorFolder + '/fpVec_' + str(tfailSec) + '.dat'

    # Make the envelope (should be low memory, most points will be thrown out as redundant)
    print 'StoreFootprintAsVector'
    myFootprint = PyFootprint(curSkyGrid)
    myFootprint.StoreFootprintAsVector(outfileStr)

    # # curSkyGrid.DumpGridToMatlab(footprintVectorFolder + '/matlabProbability_' + str(int(tfailSec)) + '.m')

    # Return the goods
    return EV_strike, outfileStr




def HaversineDistance(originDeg, destinationDeg):
    lat1, lon1 = originDeg
    lat2, lon2 = destinationDeg
    radius = 6371 # km
    
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    
    return d


### CURRENTLY IN USE (Called by falcon9.py)
def PlotNominalTrajectories(profiles, curMission, limitSec):    
    # Plot the trajectory in Google Earth
    import orbitTools as OT         ## needed for ECI2latlonalt

    # Pull out the quantities needed to plot
    stateVecStorage = profiles['stateVecStorage']
    thetagStorage = profiles['thetagStorage']

    # first index is over wind samples, second index is over actual trajectories
    numWindSamples = len(stateVecStorage)
    numTrajSamples = len(stateVecStorage[0])

    # These are storage for building up the arrays to plot
    latlonaltNumTimeSteps = []
    latlonaltStorage = np.array([],dtype=np.float64)

    # Grab some of the needed mission params
    PlanetModel = curMission['planetModel']
    deltaTsec = curMission['deltaT']
    
    # Just looking at windSample 0 for now
    for trajIX in range(numTrajSamples):
        
        # Figure out how many timesteps to plot
        numTimeSteps = 0
        numTimeStepsMax = len(thetagStorage[0][trajIX])
        
        if (limitSec < 0):
            numTimeSteps = numTimeStepsMax
        else:
            numTimeSteps = int(np.round(limitSec/deltaTsec))
            
            if (numTimeSteps > numTimeStepsMax):
                print 'ERROR!!! Youre trying to print more timesteps than there are.  Printing MAX instead'
                numTimeSteps = numTimeStepsMax
    
        latlonaltNumTimeSteps.append(numTimeSteps)
        
        yvec = np.transpose(stateVecStorage[0][trajIX][0:3])
        
        for timeIX in range(numTimeSteps):
            y = yvec[timeIX]
            thetag = thetagStorage[0][trajIX][timeIX]
            tout = 0    #if current time is zero, then can pass in thetag values at each timestep
            latlonalt = OT.ECI2latlonalt(y,thetag,PlanetModel)
            latlonaltStorage = np.concatenate((latlonaltStorage, latlonalt))
            latlonaltStorage = np.concatenate((latlonaltStorage, [0,0,0,])) # convertTJC expects six elements per line

    # Now ship this off to a google earth routine
    import data2GE
    GEfile = 'GE_Traj.kml'
    print latlonaltNumTimeSteps
    print len(latlonaltStorage)
    # data2GE.convertTJC(GEfile, latlonaltStorage, latlonaltNumTimeSteps, len(latlonaltNumTimeSteps))
    data2GE.convertTJC(GEfile, latlonaltStorage, latlonaltNumTimeSteps, 1)
    return latlonaltStorage, latlonaltNumTimeSteps



### CURRENTLY IN USE (Called by falcon9.py)
def PlotDebrisFromExplodeTime(mission, profiles, tfail, cutoffNAS = True):
    # Plot some debris trajectories that have already been generated

    # Figure out the appropriate timestep (nearestIX)
    tfailStorage = profiles['tfailStorage']
    nearestIX = (np.abs(tfailStorage[0][0]-tfail)).argmin() #Assuming first windprofile first trajectory is representative

    # open up the tProactive debris and plot it in GE
    #folder = 'testFolder'          # This will throw an error, should pass in mission1 and get debrisPickleFolder
    debrisPickleFolder = mission['debrisPickleFolder']
    GeneratedFilesFolder = mission['GeneratedFilesFolder']
    
    input = open(debrisPickleFolder + '/mpc_' + str(tfail) + '.pkl', 'rb')
    cur_mpc = pickle.load(input)
    input.close()

    # Now ship this off to a google earth routine
    from FriscoLegacy import data2GE
    GEfile = GeneratedFilesFolder + '/GE_Debris_' + str(tfail) + '.kml'
    
    deltaT = mission['deltaT']

    maxTimeSteps = np.ceil(mission['reactionTimeSeconds']/deltaT)
    data2GE.convertTJC(GEfile, cur_mpc['flatPointArray'], cur_mpc['numTimeSteps'], len(cur_mpc['numTimeSteps']), cutoffNAS, maxTimeSteps)
    # data2GE.convertTJC(GEfile, cur_mpc['flatPointArray'], cur_mpc['numTimeSteps'], len(cur_mpc['numTimeSteps']), cutoffNAS)

# ### CURRENTLY IN USE (Called by falcon9.py)
# def FindStateTimeForProactiveArchitecture(curMission, profiles, maxTime):
#     coeffIX = []     # Last index should be piece with highest ballistic coefficient, MUST BE ARRAY!!!  AAAARRRRGGGGHHHHH
# #    numTrajSamples = 20
# #    numWindSamples = 2
#     numPiecesPerSample = 5
#
#     tfail = 160             # THIS IS IN SECONDS!!!
#     shortestTime = 0;
#     shortestTimeRecord = []
#     while (shortestTime < 5.) and (tfail < maxTime):
#         tfail += 1.
#         mpc = MonteCarlo_at_tfail(curMission, coeffIX, tfail, numPiecesPerSample, profiles)
#         try:
#             shortestTime = FindShortestTimeAboveNAS(curMission, mpc)  #THIS IS IN MINUTES!!!
#         except:
#             shortestTime = -1.
#         shortestTimeRecord.append(shortestTime)
#
#         print 'tfail in seconds = ' + str(tfail) + ' gives Shortest time in minutes = ' + str(shortestTime)
#
#     return np.min(shortestTimeRecord)


### CURRENTLY IN USE (Called by falcon9.py)
def FindStateTimeForProactiveArchitecture(curMission, profiles, minTime, maxTime):
    coeffIX = []     # Last index should be piece with highest ballistic coefficient, MUST BE ARRAY!!!  AAAARRRRGGGGHHHHH
#    numTrajSamples = 20
#    numWindSamples = 2
    numPiecesPerSample = 5
    debrisPickleFolder = curMission['debrisPickleFolder']
    deltaTFail = curMission['deltaTFail']



    tfail = minTime             # THIS IS IN SECONDS!!!
    shortestTime = 0;
    shortestTimeRecord = []
    while (shortestTime < 5.5) and (tfail < maxTime):
        tfail += deltaTFail
        # mpc = MonteCarlo_at_tfail(curMission, coeffIX, tfail, numPiecesPerSample, profiles)

	if curMission['isReentry'] == True:
	    # Hack, just looking at first wind profile for now.
	    inFileName = debrisPickleFolder + '/mpc_0_' + str(int(tfail)) + '.pkl'
	else:
	    inFileName = debrisPickleFolder + '/mpc_' + str(tfail) + '.pkl'

#        input = open(debrisPickleFolder + '/mpc_' + str(tfail) + '.pkl', 'rb')
        input = open(inFileName)
        mpc = pickle.load(input)
        input.close()

        try:
            shortestTime = FindShortestTimeAboveNAS(curMission, mpc)  #THIS IS IN MINUTES!!!
        except:
            shortestTime = -1.
        shortestTimeRecord.append(shortestTime)

        print 'tfail in seconds = ' + str(tfail) + ' gives Shortest time in minutes = ' + str(shortestTime)
        print 'Maximum time of debris = {0}\n'.format(max(mpc['numTimeSteps']/60.))

    return np.max(shortestTimeRecord)


### CURRENTLY IN USE (Called by falcon9.py --> FindStateTimeForProactiveArchitecture)
def FindShortestTimeAboveNAS(curMission, mpc):
    numTimeSteps = mpc['numTimeSteps']
    flatPointArray = mpc['flatPointArray']
    big = mpc['sizeFlatPointArray']
    dt = curMission['deltaT']

    if (len(flatPointArray) == 0):
        raise IndexError

    topOfNAS = 18288    #meters

    shortestTime = 1e7  #minutes
    prevNumStepsSum = 0
    for tx in range(0,numTimeSteps.size):
        
        # The number of steps for the trajectory tx
        curNumSteps = numTimeSteps[tx]
        
        # Extract and reshape the vector of altitudes for this trajectory
        curAltitude = flatPointArray[big*prevNumStepsSum:big*(prevNumStepsSum + curNumSteps)].reshape(curNumSteps,big)[:,2]

        # Find the last point where the debris was above the NAS and convert it to a time
        # print 'curAltitude = ' + str(curAltitude)
        # print 'what it do0d?'
        # print 'topOfNAS = ' + str(topOfNAS)
        try:
            lastAltAboveNAS = curAltitude[curAltitude > topOfNAS][-1]
            lastTimeIndexAboveNAS = np.where(curAltitude == lastAltAboveNAS)[0][0]
            timeAboveNAS = lastTimeIndexAboveNAS * dt/60.
        except:
            timeAboveNAS = 0.

        # Save the time if it's the shortest we've seen thus far
        if (timeAboveNAS < shortestTime):
            shortestTime = timeAboveNAS

        # Increment so that we can pull the next altitude array out of the flatPointsArray
        prevNumStepsSum += curNumSteps


    return shortestTime





### CURRENTLY IN USE (Called by falcon9.py)
def GenerateWindTrajProfiles(curMission, numTrajSamples, numWindSamples):
    
    # To read in
    import pickle
    [altitudeList,densityMeanList,uMeanList,vMeanList, \
        wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist] \
        = pickle.load(open(curMission['atmospherePickle'],'rb'))
    
    # dt = curMission['deltaT']
    #
    # # These are some values that we're just going to hardcode
    # tLaunchDesired = 0.         # Delay in seconds of launch from nominal time
    # Toffset_std = 0.01;         # Standard Dev of the thrust offset coeff (mean of 1)
    
    # Get the mission list.  Eventually this will change format
#    missionList = loadmat(curMission['trajMatFile'])['missionList'] # loading trajectories  from SPOT
#    missionList = curMission

    # if len(curMission['propagationParamFile']) > 0:
    #     # We were given a thrust file to propagate
    #     doPropagation = True
    # else:
    #     # We were given a prepropagated trajectory and computed statevecs from that
    #     doPropagation = False

    # Generate all the atmosphere profiles and trajectories up front
    atmStorage = []
    stateVecStorage = []
    thetagStorage = []    #This is a constant for a given time of launch and failure
    tfailStorage = []
    
    # Generate the wind profiles first
    for windIX in range(numWindSamples):
        print 'windIX = ' + str(windIX)
        # generating new atmos profile
        densityList,uList,vList,wList = AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
        atmStorage.append([altitudeList,densityList,uList,vList,wList])     #altitude is in meters, density kg/m^3, 
        
        # stateVecStorage_TEMP = []
        # thetagStorage_TEMP = []
        # tfailStorage_TEMP = []
        # For every wind profile, generate a trajectory

        stateVecStorage_TEMP, thetagStorage_TEMP, tfailStorage_TEMP = \
            GenerateStateVectorProfilesForWindIX(curMission, numTrajSamples, atmStorage[windIX])

    #     if (doPropagation):
    #         for trajIX in range(numTrajSamples):
    #             if (numTrajSamples == 1):
    #                 # Just run the nominal case
    #                 ThrustOffsetAngDeg = [0,0] # offset in thrust angles
    #                 Toffset = 1. # value to multiply the thrust...
    #             else:
    #                 ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
    #                 Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
    #
    # #            StateVecProfile,thetag,tfail_real = GenerateStateVectorProfile(missionList, tLaunchDesired, dt, atmStorage[windIX], Toffset, ThrustOffsetAngDeg)
    #
    #             # Doesn't use dt!!!!  Probably comes back with dtval
    #             StateVecProfile,thetag,tfail_real = GenerateStateVectorProfile(curMission, dt, atmStorage[windIX], Toffset, ThrustOffsetAngDeg)
    #
    #             stateVecStorage_TEMP.append(StateVecProfile)
    #             thetagStorage_TEMP.append(thetag)
    #             tfailStorage_TEMP.append(tfail_real)
    #     else:
    #         # Just repeat the calculated nominal state vecs.  Don't allow for multiple trajectories
    #         sv = curMission['stateVecProfile']
    #         # stateVecProfile.append((curTime, thetag, x_eci, y_eci, z_eci, Vx, Vy, Vz))
    #         # StateVecProfile = sv[:,2:]
    #
    #         x_eci = sv[:,2]
    #         y_eci = sv[:,3]
    #         z_eci = sv[:,4]
    #
    #         Vx = sv[:,5]
    #         Vy = sv[:,6]
    #         Vz = sv[:,7]
    #
    #         VRelMag = sv[:,8]
    #
    #         thetag = sv[:,1]
    #         tfail_real = sv[:,0]
    #         stateVecStorage_TEMP.append([x_eci, y_eci, z_eci, Vx, Vy, Vz, VRelMag])
    #         thetagStorage_TEMP.append(thetag)
    #         tfailStorage_TEMP.append(tfail_real)

        
        stateVecStorage.append(stateVecStorage_TEMP)
        thetagStorage.append(thetagStorage_TEMP)
        tfailStorage.append(tfailStorage_TEMP)
    
    return atmStorage, stateVecStorage, thetagStorage, tfailStorage









### CURRENTLY IN USE (Called by Columbia.py)
def GenerateWindTrajProfilesDirectional(curMission, numTrajSamples, numWindSamples, angleLow, angleHi, windMagCoeff):
    """Appears to take in an atmoprofile, disregards the standardDev values, and instead uniformly distributes the directional variation uniformly between angleLow and Hi."""
    # To read in
    import pickle
    [altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist] \
        = pickle.load(open(curMission['atmospherePickle'],'rb'))

    # Generate all the atmosphere profiles and trajectories up front
    atmStorage = []
    stateVecStorage = []
    thetagStorage = []    #This is a constant for a given time of launch and failure
    tfailStorage = []

    # print np.where(altitudeList < 11)[0]
    # print altitudeList[np.where(altitudeList < 11)[0]]
    # # print zip(altitudeList[:,0],densityMeanList[:,0])
    # for ix in range(len(altitudeList)):
    #     print "{0}, {1}".format(altitudeList[ix], densityMeanList[ix])

    # This index is around 10m
    around10mIX = np.where(altitudeList < 11)[0][0]

    # Generate the wind profiles first
    for windIX in range(numWindSamples):
        print 'windIX = ' + str(windIX)

        # ==== Normalizing the wind profile about 10m alt ====
        # Convert to np array for coming math
        uListCol = np.array(uMeanList)
        vListCol = np.array(vMeanList)
        wListCol = np.array(wMeanList)

        # Get the profile of magntitudes all the way up
        nominalWindMagProfile = np.sqrt(uListCol**2 + vListCol**2 + wListCol**2)


        #
        # Determine the direction of the wind
        betaWind    = np.random.uniform(angleLow,angleHi)
        #

        # # Get the magnitude at this index
        # normVal = nominalWindMagProfile[around10mIX][0]
        # # Normalize the wind magnitude profile referenced to about 10m, then multiply by the velocity supplied
        # VwindMag = (nominalWindMagProfile/normVal)*windMag

        VwindMag = nominalWindMagProfile * windMagCoeff

        # Parse the wind into directions
        VeastWind   = VwindMag*np.cos(betaWind*np.pi/180.)
        VnorthWind  = VwindMag*np.sin(betaWind*np.pi/180.)

        # Load the wind into the lists that will get stored
        uList       = VeastWind
        vList       = VnorthWind
        wList       = wListCol      #np.max(wMeanList) = 0.39029202606081859, but is zero mostly.  Just add this on.


        # # # zeroing out the wind profile
        # # uListCol = 0.0*np.array(uMeanList)
        # # vListCol = 0.0*np.array(vMeanList)
        # # wListCol = 0.0*np.array(wMeanList)
        #
        # # setting wind direction on top of it
        # betaWind    = np.random.uniform(angleLow,angleHi)
        # VwindMag    = windMag
        # VeastWind   = VwindMag*np.cos(betaWind*np.pi/180.)
        # VnorthWind  = VwindMag*np.sin(betaWind*np.pi/180.)
        # uList       = uListCol + VeastWind
        # vList       = vListCol + VnorthWind
        # wList       = wListCol

        # and just using the mean density because that's what the initial state vector is tuned to
        densityList = densityMeanList

        # atmStorage = []
        atmStorage.append([altitudeList,densityList,uList,vList,wList])     #altitude is in meters, density kg/m^3,

        # # generating new atmos profile
        # densityList,uList,vList,wList = AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
        # atmStorage.append([altitudeList,densityList,uList,vList,wList])     #altitude is in meters, density kg/m^3,

        # Now get the trajectories
        stateVecStorage_TEMP, thetagStorage_TEMP, tfailStorage_TEMP = \
            GenerateStateVectorProfilesForWindIX(curMission, numTrajSamples, atmStorage[windIX])

        stateVecStorage.append(stateVecStorage_TEMP)
        thetagStorage.append(thetagStorage_TEMP)
        tfailStorage.append(tfailStorage_TEMP)

    return atmStorage, stateVecStorage, thetagStorage, tfailStorage











def GenerateStateVectorProfilesForWindIX(curMission, numTrajSamples, atmProfile):

    # These are some values that we're just going to hardcode
    tLaunchDesired = 0.         # Delay in seconds of launch from nominal time
    Toffset_std = 0.01;         # Standard Dev of the thrust offset coeff (mean of 1)
    dt = curMission['deltaT']
    planetModel = curMission['planetModel']

    if len(curMission['propagationParamFile']) > 0:
        # We were given a thrust file to propagate
        doPropagation = True
    else:
        # We were given a prepropagated trajectory and computed statevecs from that
        doPropagation = False

    stateVecStorage_TEMP = []
    thetagStorage_TEMP = []
    tfailStorage_TEMP = []

    if (doPropagation):
        # print 'YOU NEED TO UPDATE THIS FUNCTION TO INCLUDE cloption and loverd options.  EXITING!!!'
        # sys.exit()
        for trajIX in range(numTrajSamples):
            if (numTrajSamples == 1):
                # Just run the nominal case
                ThrustOffsetAngDeg = [0,0] # offset in thrust angles
                Toffset = 1. # value to multiply the thrust...
            else:
                ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
                Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...

            # Doesn't use dt!!!!  Probably comes back with dtval
            StateVecProfile,thetag,tfail_real = GenerateSingleStateVectorProfile(curMission, dt, atmProfile, Toffset, ThrustOffsetAngDeg)

            stateVecStorage_TEMP.append(StateVecProfile)
            thetagStorage_TEMP.append(thetag)
            tfailStorage_TEMP.append(tfail_real)
    else:
        # Just repeat the calculated nominal state vecs.  Don't allow for multiple trajectories
        sv = curMission['stateVecProfile']
        # stateVecProfile.append((curTime, thetag, x_eci, y_eci, z_eci, Vx, Vy, Vz))
        # StateVecProfile = sv[:,2:]

        x_eci = sv[:,2]
        y_eci = sv[:,3]
        z_eci = sv[:,4]

        Vx = sv[:,5]
        Vy = sv[:,6]
        Vz = sv[:,7]

        VRelMag = sv[:,8]

        thetag = sv[:,1]
        tfail_real = sv[:,0]
        stateVecStorage_TEMP.append([x_eci, y_eci, z_eci, Vx, Vy, Vz, VRelMag])
        thetagStorage_TEMP.append(thetag)
        tfailStorage_TEMP.append(tfail_real)

        # if numTrajSamples == 1, this will return [] and do nothing
        if numTrajSamples > 1:
            import orbitTools as ot
            print "NOT YET IMPLEMENTED!!!\n"
            sys.exit()

            # Convert state vectors back to lat lon alts
            for cursv in sv:
                thetag = cursv[1]
                r_eci = cursv[2:5]
                v_eci = cursv[5:8]

                ot.ECI2latlonalt(r_eci, thetag, planetModel)



            for trajIX in range(1,numTrajSamples):
                # Okay, we were asked for multiple trajectories, let's put some randomness in there
                None




    return stateVecStorage_TEMP, thetagStorage_TEMP, tfailStorage_TEMP



### CURRENTLY IN USE (Called by falcon9.py --> GenerateWindTrajProfiles)
def GenerateSingleStateVectorProfile(missionList, dt, atmProfile, Toffset, ThrustOffsetAngDeg):
    
#    myTrajList = trajWHATHWHATWHAT.generateRandomTraj(missionList,tLaunchDesired,dt,atmProfile,Toffset,ThrustOffsetAngDeg)
    
    myTrajList = getStateVector(missionList, atmProfile, Toffset, ThrustOffsetAngDeg)   # I think this comes back with timesteps of dtval
    
    numStages = len(myTrajList)
    
    # Should have made a note of this when I did it the first time...but I think here dstack is switching when you index by the stage.
    #   For example, the xvals for the first stage are myTrajList[0][0] and for the second stage are myTrajList[1][0].  But after dstack,
    #   you'll get a vector like myTrajStack[0][1] = (myTrajList[0][0], myTrajList[1][0])...so this effectively moves the stage index from
    #   position 1 to position 2
    myTrajStack = np.dstack(myTrajList)
    myTrajStackFinal = []
    
    if (numStages == 2):
        for ix in range(9):
            myTrajStackFinal.append(np.hstack( [myTrajStack[0][ix][0], myTrajStack[0][ix][1][1:]]  ))
    elif (numStages == 3):
        for ix in range(9):
            myTrajStackFinal.append(np.hstack( [myTrajStack[0][ix][0], myTrajStack[0][ix][1][1:], myTrajStack[0][ix][2][1:]]  ))
    else:
        print 'ERROR!  GenerateSingleStateVectorProfile needs to be programmed to handle the case of stages = ' + str(numStages)

    [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector] = myTrajStackFinal    #index 0=first stage, 1=coast period, 2=second stage
    stateVector = [xc,yc,zc,Vxc,Vyc,Vzc]
    
    return stateVector, thetagVector, tc


#### CURRENTLY IN USE (Called by falcon9.py)
#def MonteCarlo_until_tfail(curMission, profiles, tfail, folder):
#    
#    # Make sure that the output directory exists
#    folderPath = os.path.abspath(folder)
#    if not os.path.exists(folderPath):
#        os.makedirs(folderPath)
#    
#    coeffIX = []     # Empty array means use all ballistic coeffs
#    numPiecesPerSample = 5
#    
#    curTime = tfail*1.0            # THIS IS IN SECONDS!!!
#    while (curTime > 0):
#        MonteCarlo_at_tfail_and_record(curMission, coeffIX, curTime, numPiecesPerSample, profiles, folderPath)
#        curTime -= 5.

# ### CURRENTLY IN USE (Called by falcon9.py)
# def MonteCarlo_until_tfail(curMission, profiles, tfail):
#
#     coeffIX = []     # Empty array means use all ballistic coeffs
#     numPiecesPerSample = curMission['numPiecesPerSample']
#     deltaTFail = curMission['deltaTFail']
#
#     # Do the parallel
#     import pp
#     job_server = pp.Server()
#     print 'number of cpus = ' + str(job_server.get_ncpus())
#
#     # by placing -deltaTFail as the lower limit, we will include time = zero
#     timeVec = np.arange(tfail*1.0,-deltaTFail,-deltaTFail)        #curTime is in seconds
#     # timeVec = np.arange(tfail*1.0,0,-deltaTFail)        #curTime is in seconds
#     jobs = [job_server.submit(MonteCarlo_at_tfail_and_record, \
#                               args=(curMission, coeffIX, curTime, numPiecesPerSample, profiles), \
#                               depfuncs=(MonteCarlo_at_tfail,), \
#                               modules=('numpy as np','debrisReader as DR', 'orbitTools', 'debrisPropagation as dp'), \
#                               callback=finished001) for curTime in timeVec]
#
#     job_server.wait()
#     job_server.print_stats()
#     for job in jobs:
#         print job()


def MonteCarlo_until_tfail(curMission, profiles, tlow, thigh):

    coeffIX = []     # Empty array means use all ballistic coeffs
    numPiecesPerSample = curMission['numPiecesPerSample']
    deltaTFail = curMission['deltaTFail']

    # by placing -deltaTFail as the lower limit, we will include time = zero
    timeVec = np.arange(thigh*1.0,tlow-deltaTFail,-deltaTFail)        #curTime is in seconds

    # Test one
    print 'TESTING'
    MonteCarlo_at_tfail_and_record(curMission, coeffIX, timeVec[-1], numPiecesPerSample, profiles)
    print 'DONE TESTING'

    # Do the parallel
    # Do the parallel
    jobs = []   # Holds the results of the envelope creation
    if curMission['numNodes'] == 1:
        # This means we're running into memory or malloc issues with pp, so don't use it
        # jobs = [genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) for ix in range(len(timeRange))]
        jobs = [MonteCarlo_at_tfail_and_record(curMission, coeffIX, curTime, numPiecesPerSample, profiles) \
                for curTime in timeVec]
    else:
        import pp
        job_server = pp.Server()
        job_server.set_ncpus(curMission['numNodes'])
        print 'number of cpus = ' + str(job_server.get_ncpus())

        jobs = [job_server.submit(MonteCarlo_at_tfail_and_record, \
                                  args=(curMission, coeffIX, curTime, numPiecesPerSample, profiles), \
                                  depfuncs=(MonteCarlo_at_tfail,), \
                                  modules=('numpy as np','from FriscoLegacy import debrisReader as DR', 'from FriscoLegacy import orbitTools', 'from FriscoLegacy import debrisPropagation as dp'), \
                                  callback=finished001) for curTime in timeVec]

        print "This will fail because the modules like debrisReader are now inside a package.  \
                Must include package info or put all of these modules in the root package. 1182 "
        # sys.exit()

        job_server.wait()
        job_server.print_stats()
        for job in jobs:
            print job()

def finished001(val):
    print "done with " + str(val)


### Wrapper function used for efficient parallelization of MonteCarlo_until_tfail (Called by MonteCarlo_until_tfail)
def MonteCarlo_at_tfail_and_record(curMission, coeffIX, curTime, numPiecesPerSample, profiles):
    #print 'MonteCarloing at tstep = ' + str(curTime)
    print 'going in'
    mpc = MonteCarlo_at_tfail(curMission, coeffIX, curTime, numPiecesPerSample, profiles)
    print 'made it out'
    debrisPickleFolder = curMission['debrisPickleFolder']

    # print 'COMMENTED OUT WRITING TO FILE FOR DEBUGGING PURPOSES'
    # Make sure that the output directory exists
    folderPath = os.path.abspath(debrisPickleFolder)
    if not os.path.exists(folderPath):
        os.makedirs(folderPath)

    output = open(folderPath + '/mpc_' + str(curTime) + '.pkl', 'wb')
    pickle.dump(mpc,output,2)
    output.close()

    return curTime

### CURRENTLY IN USE (Called by falcon9.py --> MonteCarlo_until_tfail AND FindStateTimeForProactiveArchitecture)
def MonteCarlo_at_tfail(curMission, coeffIX, tfail, numPiecesPerSample, profiles):
    # Unpack some of the mission variables
    debrisCatalogFile = curMission['debrisCatPath'] + curMission['debrisCatFile']
    debrisCatalogFilePATH = curMission['debrisCatPath']
    dt = curMission['deltaT']
    ndtinterval = curMission['debrisTimeLimitSec']/dt
    
    # ndtinterval = 1*3600/dt     # is the upper bound for how long you think the debris propagation will run (5hrs)
    tLaunchDesired = 0.         # Delay in seconds of launch from nominal time
    
    # Unpack the profiles
    atmStorage      = profiles['atmStorage']
    thetagStorage   = profiles['thetagStorage']
    stateVecStorage = profiles['stateVecStorage']
    tfailStorage    = profiles['tfailStorage']
    
    # Find sizes of profiles
    numWindSamples  = len(atmStorage)
    numTrajSamples  = len(stateVecStorage[0])
    
    # Need Vmag to find correct debris catalog
    if (tfail > tfailStorage[0][0][-1]):
        # We don't have information out this far.  Kill it
        mpc = dict(flatPointArray = [])
        return mpc

    nearestIX = (np.abs(tfailStorage[0][0]-tfail)).argmin() #Assuming first windprofile first trajectory is representative

    x = stateVecStorage[0][0][0][nearestIX]
    y = stateVecStorage[0][0][1][nearestIX]
    z = stateVecStorage[0][0][2][nearestIX]
    Vx = stateVecStorage[0][0][3][nearestIX]
    Vy = stateVecStorage[0][0][4][nearestIX]
    Vz = stateVecStorage[0][0][5][nearestIX]
    
    omegaE = curMission['omegaE']
    Vrot = np.cross([0,0,omegaE],[x,y,z])
    Vinf = np.array([Vx,Vy,Vz]) - Vrot
    Vmag = (Vinf[0]**2 + Vinf[1]**2 + Vinf[2]**2)**.5
    
    # getting debris catalog information, getting desired debris catalog group
    velcat,timecat,catalogList = DR.readDebris(debrisCatalogFile)
    desiredCatalog = DR.desiredList(velcat,timecat,catalogList,Vmag,tfail)
    
    # generating random debris pieces
    massList,arefList,velList,ncdList,minfcdList,CDList,nclList,minfclList,CLList,arefMean,numberOfPiecesMeanList,blast=DR.generateDebris(desiredCatalog,numPiecesPerSample,debrisCatalogFilePATH)
    
    # Only want to do a single debris coefficient, use the one passed in
    if (coeffIX == []):
        # If it's empty, then let's look at all the debris pieces
        coeffIX = range(len(massList))
    # If it's not empty, then just look at the pieces specified (probably -1)
    
    # These two arrays will be appended to and returned as the solution of this function
    debrisNumTimeSteps = []
    # debrisStorage = np.array([],dtype=np.float64)
    debrisID = []
    debrisMass = []
    debrisArea = []
    testStorage = []
    
    for windIX in range(numWindSamples):
        altitudeList = atmStorage[windIX][0]
        densityList = atmStorage[windIX][1]
        uList = atmStorage[windIX][2]
        vList = atmStorage[windIX][3]
        wList = atmStorage[windIX][4]
        
        for trajIX in range(numTrajSamples):
            #            stateVec = stateVecStorage[windIX][trajIX]
            thetag = thetagStorage[windIX][trajIX][nearestIX]
            
            x = stateVecStorage[windIX][trajIX][0][nearestIX]
            y = stateVecStorage[windIX][trajIX][1][nearestIX]
            z = stateVecStorage[windIX][trajIX][2][nearestIX]
            Vx = stateVecStorage[windIX][trajIX][3][nearestIX]
            Vy = stateVecStorage[windIX][trajIX][4][nearestIX]
            Vz = stateVecStorage[windIX][trajIX][5][nearestIX]
            stateVec = [x,y,z,Vx,Vy,Vz]
            
            for debrisIndex in coeffIX:
                ncdval = ncdList[debrisIndex]
                minfcdval = minfcdList[debrisIndex]
                nclval = nclList[debrisIndex]
                minfclval = minfclList[debrisIndex]
                CDval = CDList[debrisIndex]
                CLval = CLList[debrisIndex]

                # print '[wind][traj][beta] = ' + '[{0}][{1}][{2}]'.format(windIX, trajIX, debrisIndex)
                
                for pieceIndex in range(numPiecesPerSample):
                    
                    massval = massList[debrisIndex][pieceIndex]
                    arefval = arefList[debrisIndex][pieceIndex]
                    velImpval = velList[debrisIndex][pieceIndex]
                    #velImpval = 100.  # addition...just to check behavior
                    theta1 = np.random.uniform(0.0,2*np.pi)
                    theta2 = np.random.uniform(0.0,2*np.pi) #it is not required this to be 0 to 2pi, domain is covered with 0 to pi since prev angle helps cover the entire sphere
                    Vdeb = np.array([[velImpval,0,0]]).T # velocity impulse calculation
                    Rz = orbitTools.cRnMatrix_Zaxis(theta1)
                    Rx = orbitTools.cRnMatrix_Xaxis(theta2)
                    Vdeb = np.dot(Rz,Vdeb)
                    Vdeb = np.dot(Rx,Vdeb)
                    Vdeb = Vdeb[:,0] # making it a 1 by 3 vector...impulse velocity in cartesin coordinates
                    #print arefval
                    # debris propagation calculation
                    
                    #            print 'pieceIX = ' + str(pieceIndex) + '   mass = ' + str(massval)
                                        
                    debrisResults, numFinalSteps = dp.debrispropagation(initialstate = stateVec,
                                                                        debrisvel = Vdeb,
                                                                        mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=1,
                                                                        minfcl=minfclval,cl=CLval,loverd = 0,atmosoption=2,altitudelist=altitudeList,
                                                                        densitylist=densityList,ulist=uList,vlist=vList,
                                                                        wlist=wList,geoptions=0,filename='none',planetmodel=0,dtinterval = dt,ndtinterval = ndtinterval,thetag0=thetag,
                                                                        ncd=ncdval,
                                                                        ncl=nclval,
                                                                        nlist=len(altitudeList))                    
                    #                pdb.set_trace()
                    
#                    print 'QUITTING AFTER ONE ITERATION AND RETURNING TEST VALUE FOR DEBUGGING PURPOSES!!!!'
#                    return debrisResults, numFinalSteps
                    
                    # Explanation of these manipulations:
                    # debrisResults is returned as a fixed-length one-dimensional array.  Fortran makes it way larger than it
                    #   would ever need to be because it doesn't have the ability to do dynamic memory allocation and it doesn't
                    #   know at the beginning of the propagation how long it will take the debris to hit the ground.  Thus, we need
                    #   to trim all the unused stuff off of the end
                    #   * debrisResults[0:(numFinalSteps)*sizeResults]
                    # Then we want to unflatten the array into a multidimensional array
                    #   * .reshape(numFinalSteps,sizeResults)
                    # But there might be return values that we don't care about, so select only the columns that we want to store.
                    #   In my case, that's the first six values (position and velocity)
                    #   * [:,:6]
                    # Later when we store the values in debrisStorage, we'll flatten them out and make this is a massive 1D array
                    #   * .flatten()
                    sizeResults = 7
                    sizeKeeping = 6
                    debrisResults = debrisResults[0:(numFinalSteps)*sizeResults].reshape(numFinalSteps,sizeResults)[:,:sizeKeeping]
                    
                    altitudeFinal = debrisResults[-1,2]
                    # if altitudeFinal>0:
                    #     print 'Warning...debris is not on the ground. Final altitude is ',altitudeFinal
                                        
                    # Saves the stuff I care about
                    # For some reason, those extra parentheses are necessary to avoid errors when concatenating the first (empty) array
                    # debrisStorage = np.concatenate(( debrisStorage,  debrisResults.flatten() ))
                    testStorage.append(debrisResults.flatten() )
                    debrisNumTimeSteps.append(numFinalSteps)
                    debrisID.append(debrisIndex)
                    debrisMass.append(massval)
                    debrisArea.append(arefval)

    # Convert the time steps array to be numpy int32 array
    debrisNumTimeSteps = np.array(debrisNumTimeSteps,dtype=np.int32)
    debrisID = np.array(debrisID,dtype=np.int32)
    debrisMass = np.array(debrisMass,dtype=np.double)
    debrisArea = np.array(debrisArea,dtype=np.double)

    numberOfPiecesMeanList = np.array(numberOfPiecesMeanList,dtype=np.int32)

    # This is no longer needed, zero it out and eventually remove it
    arefMean = np.array(arefMean,dtype=np.double)
    arefMean = np.array([])

    # This flattens the python list testStorage and saves it as the proper numpy array
    debrisStorage = np.array([item for sublist in testStorage for item in sublist],dtype=np.float64)


    ###### Now to read these points into a PointCloud, for input into a footprint and export to GE
    UTC = curMission['initialUTC'] + (tLaunchDesired + tfail)/(24*3600.)
    # zBinHeightKm = 5.0

    numRuns = len(debrisNumTimeSteps)       # This is the total number of pieces of debris simulated here
    maxTime = debrisNumTimeSteps.max()*dt
    
    # Pack things up in a dictionary (for pickling)
    mpc = dict(flatPointArray = debrisStorage, numPieces = numRuns, numTimeSteps = debrisNumTimeSteps, maxTime = maxTime, deltaTsec = curMission['deltaT'],
               UTC = UTC, all_points_delta_t = curMission['all_points_delta_t'],
               launchLat = curMission['launchLat'], launchLon = curMission['launchLon'], launchAzimuth = curMission['launchAzimuth'],
               debrisID = debrisID, arefMeanList = arefMean, numberOfPiecesMeanList = numberOfPiecesMeanList,
               debrisMass = debrisMass, debrisArea = debrisArea, sizeFlatPointArray = sizeKeeping)
            
    return mpc






def MonteCarlo_Distributed_Reentry_Wrapper_CAIB(curMission, coeffIX, numPiecesPerSample, lowerTime, upperTime, profiles):

    # Unpack some of the mission variables
    debrisCatalogFile = curMission['debrisCatPath'] + curMission['debrisCatFile']
    # Need Vmag to find correct debris catalog
    Vmag    = 1    #No dependence on this within the catalog

    # tfail   = 1     # time is now very important
    # # getting debris catalog information, getting desired debris catalog group
    velcat,timecat,catalogList = DR.readDebris(debrisCatalogFile)
    catalogList = zip(timecat,catalogList)  #packaging these two together

    if (timecat[-1] > profiles['tfailStorage'][0][0][-1]):
        # We don't have information out this far.  Kill it
        mpc = dict(flatPointArray = [])
        print 'Returning an empty mpc'
        print '{0} > {1}'.format(timecat[-1], profiles['tfailStorage'][0][0][-1])
        return mpc

    # desiredCatalog = DR.desiredList(velcat,timecat,catalogList,Vmag,tfail)

    # # Only want to do a single debris coefficient, use the one passed in
    # if (coeffIX == []):
    #     # If it's empty, then let's look at all the debris pieces
    #     coeffIX = range(len(desiredCatalog))
    # # If it's not empty, then just look at the pieces specified (probably -1)

    # print 'coeffIX = ' + str(coeffIX)

    numWindSamples = len(profiles['atmStorage'])

    # Use the passed-in values
    lowerBreakLimit = lowerTime
    upperBreakLimit = upperTime
    # # Unless the time range is not specified, then explode everywhere
    # if (len(lowerTime) == 0) and (len(upperTime) == 0):
    #     # Set the times to explode to coincide with the nominal profile
    #     lowerBreakLimit = profiles['tfailStorage'][0][0][0] #Assuming first trajectory is representative
    #     upperBreakLimit = profiles['tfailStorage'][0][0][-1]


    deltaTFail = curMission['deltaTFail']
    # Won't include very last time step which may well be negative altitude and thus undesireable
    timeVec = np.arange(lowerBreakLimit,upperBreakLimit,deltaTFail)        #curTime is in seconds
    # tfail = 300


    # ## Commenting out for now, but will be resurrected after Columbia stuff is sorted
    # # for windIX in range(numWindSamples):
    # windIX = 0
    # cur_mpc = MonteCarlo_Distributed_Reentry_and_Record_CAIB(curMission, catalogList, coeffIX, numPiecesPerSample,
    #                                          profiles, lowerBreakLimit, upperBreakLimit, windIX)


    # Do the parallel
    jobs = []   # Holds the results of the envelope creation
    if curMission['numNodes'] == 1:
        # This means we're running into memory or malloc issues with pp, so don't use it
        # jobs = [genFootprint(mission1, timeRange[ix], pFailThisTimestepVec[ix]) for ix in range(len(timeRange))]
        jobs = [\
                MonteCarlo_Distributed_Reentry_and_Record_CAIB(curMission, catalogList, coeffIX, numPiecesPerSample, \
                                                     profiles, lowerBreakLimit, upperBreakLimit, tfail, windIX) \
                for windIX in range(numWindSamples) for tfail in timeVec]
    else:
        import pp
        job_server = pp.Server()
        job_server.set_ncpus(curMission['numNodes'])
        print 'number of cpus = ' + str(job_server.get_ncpus())

        # by placing -deltaTFail as the lower limit, we will include time = zero
        # timeVec = np.arange(thigh*1.0,tlow-deltaTFail,-deltaTFail)        #curTime is in seconds
        # timeVec = np.arange(tfail*1.0,-deltaTFail,-deltaTFail)        #curTime is in seconds
        # timeVec = np.arange(tfail*1.0,0,-deltaTFail)        #curTime is in seconds

        # I think maybe this is just too big.  Split up the parallelization
        for tfail in timeVec:
            jobs = [job_server.submit(MonteCarlo_Distributed_Reentry_and_Record_CAIB, \
                                      args=(curMission, catalogList, coeffIX, numPiecesPerSample, \
                                                         profiles, lowerBreakLimit, upperBreakLimit, tfail, windIX), \
                                      depfuncs=(MonteCarlo_Distributed_Reentry_CAIB,), \
                                      modules=('numpy as np','from FriscoLegacy import debrisReader as DR', 'from FriscoLegacy import orbitTools', 'from FriscoLegacy import debrisPropagation as dp'), \
                                      callback=finishedDistributed) for windIX in range(numWindSamples)]  
            job_server.wait()
            job_server.print_stats()
            for job in jobs:
                print job()
        job_server.destroy()


        # jobs = [job_server.submit(MonteCarlo_Distributed_Reentry_and_Record_CAIB, \
        #                           args=(curMission, catalogList, coeffIX, numPiecesPerSample, \
        #                                              profiles, lowerBreakLimit, upperBreakLimit, tfail, windIX), \
        #                           depfuncs=(MonteCarlo_Distributed_Reentry_CAIB,), \
        #                           modules=('numpy as np','from FriscoLegacy import debrisReader as DR', 'from FriscoLegacy import orbitTools', 'from FriscoLegacy import debrisPropagation as dp'), \
        #                           callback=finishedDistributed) for windIX in range(numWindSamples) for tfail in timeVec]

        # job_server.wait()
        # job_server.print_stats()
        # for job in jobs:
        #     print job()
        # job_server.destroy()

def finishedDistributed(val):
    print 'done: windIX = {0}, tfail = {1}'.format(val[0], val[1])


### Wrapper function used for efficient parallelization of MonteCarlo_until_tfail (Called by MonteCarlo_until_tfail)
def MonteCarlo_Distributed_Reentry_and_Record_CAIB(curMission, catalogList, coeffIX, numPiecesPerSample,
                                                 profiles, lowerBreakLimit, upperBreakLimit, tfail, windIX):

    # mpc = MonteCarlo_at_tfail(curMission, coeffIX, curTime, numPiecesPerSample, profiles)

    cur_mpc = MonteCarlo_Distributed_Reentry_CAIB(curMission, catalogList, coeffIX, numPiecesPerSample,
                                                 profiles, lowerBreakLimit, upperBreakLimit, tfail, windIX)

    # # print 'COMMENTED OUT WRITING TO FILE FOR DEBUGGING PURPOSES'
    # # Make sure that the output directory exists
    # folderPath = os.path.abspath(debrisPickleFolder)
    # if not os.path.exists(folderPath):
    #     os.makedirs(folderPath)

    outFileName = '{0}/mpc_{1}_{2}.pkl'.format(curMission['debrisPickleFolder'], windIX, int(tfail))
    output = open(outFileName, 'wb')
    # output = open(curMission['debrisPickleFolder'] + '/mpc_' + str(windIX) + '.pkl', 'wb')
    pickle.dump(cur_mpc,output,2)
    output.close()

    return (windIX, tfail)

    # return curTime


### CURRENTLY IN USE (Called by Columbia.py -->
def MonteCarlo_Distributed_Reentry_CAIB(curMission, catalogList, coeffIX, numPiecesPerSample, profiles, lowerBreakLimit, upperBreakLimit, tfail, windIX):

    debrisCatalogFilePATH = curMission['debrisCatPath']
    dt = curMission['deltaT']
    ndtinterval = curMission['debrisTimeLimitSec']/dt

    cloption = 1
    if curMission['useLoverD'] == True:
        cloption = 0
    loverd = curMission['loverd']

    # ndtinterval = 1*3600/dt     # is the upper bound for how long you think the debris propagation will run (1hrs)
    # ndtinterval = 6*60/dt     # is the upper bound for how long you think the debris propagation will run (1hrs)
    tLaunchDesired = 0.         # Delay in seconds of launch from nominal time

    # Unpack the profiles
    atmStorage = profiles['atmStorage']
    thetagStorage = profiles['thetagStorage']
    stateVecStorage = profiles['stateVecStorage']
    tfailStorage = profiles['tfailStorage']

    # Find sizes of profiles
    # numWindSamples = len(atmStorage)
    numTrajSamples = len(stateVecStorage[0])

    omegaE = curMission['omegaE']
    planetmodel = curMission['planetModel']
    # print 'planetmodel = ' + str(planetmodel)

    # These two arrays will be appended to and returned as the solution of this function
    debrisNumTimeSteps = []
    debrisID = []
    debrisMass = []
    debrisArea = []
    testStorage = []

    altitudeList    = atmStorage[windIX][0]
    densityList     = atmStorage[windIX][1]
    uList           = atmStorage[windIX][2]
    vList           = atmStorage[windIX][3]
    wList           = atmStorage[windIX][4]

    if curMission['noWind']:
        uList = 0.*uList
        vList = 0.*vList
        wList = 0.*wList

    numberOfPiecesMeanListStorage = []

    for trajIX in range(numTrajSamples):

        # debrisID must now change across times
        debrisIDcounter = 0

        # Iterate through the timesteps in the catalog and propagate those pieces
        for curTuple in catalogList:
            curTime     = curTuple[0]
            curCatalog  = curTuple[1]

            # nSamplesVec = []
            # for group in curCatalog:
            #     nSamplesVec.append(int(group[0]))
            # nSamplesVec = range(len(massList))

            # If you don't specify a numPiecesPerSample value, then use the exact values from the debris catalog
            nSamplesVec = []
            if len(numPiecesPerSample) == 0:
                for group in curCatalog:
                    nSamplesVec.append(int(group[0]))
            else:
                for group in curCatalog:
                    nSamplesVec.append(numPiecesPerSample[0])

            # print 'curTime = {0}, nSamplesVec = {1}'.format(curTime, nSamplesVec)
            # print curCatalog

            # generating random debris pieces
            massList,arefList,velList,ncdList,minfcdList,CDList,nclList,\
            minfclList,CLList,arefMean,numberOfPiecesMeanList,blast\
                = DR.generateDebris(curCatalog,nSamplesVec,debrisCatalogFilePATH)

            # Banishing coeffIX, it is no longer a good measure of anything
            for debrisIndex in range(len(nSamplesVec)):
                ncdval = ncdList[debrisIndex]
                minfcdval = minfcdList[debrisIndex]
                nclval = nclList[debrisIndex]
                minfclval = minfclList[debrisIndex]
                CDval = CDList[debrisIndex]
                CLval = CLList[debrisIndex]

                # print '[wind][traj][beta] = ' + '[{0}][{1}][{2}]'.format(windIX, trajIX, debrisIndex)

                for pieceIndex in range(nSamplesVec[debrisIndex]):

                    # sampling time during a two minute window
                    # tbreak = np.random.uniform(lowerBreakLimit,upperBreakLimit)
                    catalogTimeInterval = 5 # I know this is true for now, but should make this an input
                    tbreak = tfail + curTime + np.random.uniform(0,catalogTimeInterval)
                                #Timing notes:
                                # * tfail is time of failure
                                # * curTime is increment from tfail
                                # * Then adding randomness to time within increment

                    # This gets moved into loop for each piece within a generated debris catalog
                    nearestIX = (np.abs(tfailStorage[0][0]-tbreak)).argmin() #Assuming first windprofile first trajectory is representative

                    thetag = thetagStorage[windIX][trajIX][nearestIX]
                    x = stateVecStorage[windIX][trajIX][0][nearestIX]
                    y = stateVecStorage[windIX][trajIX][1][nearestIX]
                    z = stateVecStorage[windIX][trajIX][2][nearestIX]
                    Vx = stateVecStorage[windIX][trajIX][3][nearestIX]
                    Vy = stateVecStorage[windIX][trajIX][4][nearestIX]
                    Vz = stateVecStorage[windIX][trajIX][5][nearestIX]

                    stateVec = [x,y,z,Vx,Vy,Vz]

                    Vrot = np.cross([0,0,omegaE],[x,y,z])
                    Vinf = np.array([Vx,Vy,Vz]) - Vrot
                    Vmag = (Vinf[0]**2 + Vinf[1]**2 + Vinf[2]**2)**.5

                    massval = massList[debrisIndex][pieceIndex]
                    arefval = arefList[debrisIndex][pieceIndex]
                    velImpval = velList[debrisIndex][pieceIndex]
                    #velImpval = 100.  # addition...just to check behavior
                    theta1 = np.random.uniform(0.0,2*np.pi)
                    theta2 = np.random.uniform(0.0,2*np.pi) #it is not required this to be 0 to 2pi, domain is covered with 0 to pi since prev angle helps cover the entire sphere
                    Vdeb = np.array([[velImpval,0,0]]).T # velocity impulse calculation
                    Rz = orbitTools.cRnMatrix_Zaxis(theta1)
                    Rx = orbitTools.cRnMatrix_Xaxis(theta2)
                    Vdeb = np.dot(Rz,Vdeb)
                    Vdeb = np.dot(Rx,Vdeb)
                    Vdeb = Vdeb[:,0] # making it a 1 by 3 vector...impulse velocity in cartesin coordinates
                    #print arefval
                    # debris propagation calculation

                    # print 'pieceIndex = {0}, Vdeb/Vmag = {1}/{2}'.format(pieceIndex, (Vdeb[0]**2 + Vdeb[1]**2 + Vdeb[2]**2)**.5, Vmag)

                    #            print 'pieceIX = ' + str(pieceIndex) + '   mass = ' + str(massval)

                    # print 'tbreak = {0}, nearestIX = {1}'.format(tbreak, nearestIX)

                    debrisResults, numFinalSteps = dp.debrispropagation(initialstate = stateVec,
                                                                        debrisvel = Vdeb,
                                                                        mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=cloption,
                                                                        minfcl=minfclval,cl=CLval,loverd = loverd, atmosoption=2,altitudelist=altitudeList,
                                                                        densitylist=densityList,ulist=uList,vlist=vList,
                                                                        wlist=wList,geoptions=0,filename='none',planetmodel=planetmodel,dtinterval = dt,ndtinterval = ndtinterval,thetag0=thetag,
                                                                        ncd=ncdval,
                                                                        ncl=nclval,
                                                                        nlist=len(altitudeList))
                    #                pdb.set_trace()

        #                    print 'QUITTING AFTER ONE ITERATION AND RETURNING TEST VALUE FOR DEBUGGING PURPOSES!!!!'
        #                    return debrisResults, numFinalSteps

                    # Explanation of these manipulations:
                    # debrisResults is returned as a fixed-length one-dimensional array.  Fortran makes it way larger than it
                    #   would ever need to be because it doesn't have the ability to do dynamic memory allocation and it doesn't
                    #   know at the beginning of the propagation how long it will take the debris to hit the ground.  Thus, we need
                    #   to trim all the unused stuff off of the end
                    #   * debrisResults[0:(numFinalSteps)*sizeResults]
                    # Then we want to unflatten the array into a multidimensional array
                    #   * .reshape(numFinalSteps,sizeResults)
                    # But there might be return values that we don't care about, so select only the columns that we want to store.
                    #   In my case, that's the first six values (position and velocity)
                    #   * [:,:6]
                    # Later when we store the values in debrisStorage, we'll flatten them out and make this is a massive 1D array
                    #   * .flatten()
                    sizeResults = 7
                    sizeKeeping = 6
                    debrisResults = debrisResults[0:(numFinalSteps)*sizeResults].reshape(numFinalSteps,sizeResults)[:,:sizeKeeping]

                    altitudeFinal = debrisResults[-1,2]
                    # if altitudeFinal>0:
                    #     print 'Warning...debris is not on the ground. Final altitude is ',altitudeFinal

                    # print "numFinalSteps = {0}".format(numFinalSteps)
                    # Saves the stuff I care about
                    # For some reason, those extra parentheses are necessary to avoid errors when concatenating the first (empty) array
                    # debrisStorage = np.concatenate(( debrisStorage,  debrisResults.flatten() ))
                    testStorage.append(debrisResults.flatten() )
                    debrisNumTimeSteps.append(numFinalSteps)
                    debrisID.append(debrisIDcounter)
                    debrisMass.append(massval)
                    debrisArea.append(arefval)
                debrisIDcounter += 1
            numberOfPiecesMeanListStorage.extend(numberOfPiecesMeanList)

    # Convert the time steps array to be numpy int32 array
    # This info exists for each individual piece
    debrisNumTimeSteps = np.array(debrisNumTimeSteps,dtype=np.int32)
    debrisID = np.array(debrisID,dtype=np.int32)
    debrisMass = np.array(debrisMass,dtype=np.double)
    debrisArea = np.array(debrisArea,dtype=np.double)

    # This is overall information
    numberOfPiecesMeanList = np.array(numberOfPiecesMeanListStorage,dtype=np.int32)

    # This is no longer needed, zero out and eventually remove
    arefMean = np.array(arefMean,dtype=np.double)
    arefMean = []

    # This flattens the python list testStorage and saves it as the proper numpy array
    debrisStorage = np.array([item for sublist in testStorage for item in sublist],dtype=np.float64)


    ###### Now to read these points into a PointCloud, for input into a footprint and export to GE
    # UTC = curMission['initialUTC'] + (tLaunchDesired + tfail)/(24*3600.)
    UTC = curMission['initialUTC']
    # zBinHeightKm = 5.0

    numRuns = len(debrisNumTimeSteps)
    maxTime = debrisNumTimeSteps.max()*dt

    # Pack things up in a dictionary (for pickling)
    mpc = dict(flatPointArray = debrisStorage, numPieces = numRuns, numTimeSteps = debrisNumTimeSteps, maxTime = maxTime, deltaTsec = curMission['deltaT'],
               UTC = UTC, all_points_delta_t = curMission['all_points_delta_t'],
               launchLat = curMission['launchLat'], launchLon = curMission['launchLon'], launchAzimuth = curMission['launchAzimuth'],
               debrisID = debrisID, arefMeanList = arefMean, numberOfPiecesMeanList = numberOfPiecesMeanList,
               debrisMass = debrisMass, debrisArea = debrisArea, sizeFlatPointArray = sizeKeeping)

    return mpc








# def MonteCarlo_Distributed_Reentry_Wrapper(curMission, coeffIX, numPiecesPerSample, profiles):
#     upperBreakLimit = 135
#     lowerBreakLimit = 15
#
#     if (upperBreakLimit > profiles['tfailStorage'][0][0][-1]):
#         # We don't have information out this far.  Kill it
#         mpc = dict(flatPointArray = [])
#         return mpc
#
#     # Unpack some of the mission variables
#     debrisCatalogFile = curMission['debrisCatPath'] + curMission['debrisCatFile']
#     # Need Vmag to find correct debris catalog
#     Vmag    = 1    #No dependence on this within the catalog
#     tfail   = 1
#     # getting debris catalog information, getting desired debris catalog group
#     velcat,timecat,catalogList = DR.readDebris(debrisCatalogFile)
#     desiredCatalog = DR.desiredList(velcat,timecat,catalogList,Vmag,tfail)
#
#     # Only want to do a single debris coefficient, use the one passed in
#     if (coeffIX == []):
#         # If it's empty, then let's look at all the debris pieces
#         coeffIX = range(len(desiredCatalog))
#     # If it's not empty, then just look at the pieces specified (probably -1)
#
#     print 'coeffIX = ' + str(coeffIX)
#
#     numWindSamples = len(profiles['atmStorage'])
#
#     # ## Commenting out for now, but will be resurrected after Columbia stuff is sorted
#     # for windIX in range(numWindSamples):
#     #     cur_mpc = MonteCarlo_Distributed_Reentry_and_Record(curMission, desiredCatalog, coeffIX, numPiecesPerSample,
#     #                                              profiles, lowerBreakLimit, upperBreakLimit, windIX)
#
#
#
#     # Do the parallel
#     import pp
#     job_server = pp.Server()
#     # job_server.set_ncpus(6)
#     print 'number of cpus = ' + str(job_server.get_ncpus())
#
#     # by placing -deltaTFail as the lower limit, we will include time = zero
#     # timeVec = np.arange(thigh*1.0,tlow-deltaTFail,-deltaTFail)        #curTime is in seconds
#     # timeVec = np.arange(tfail*1.0,-deltaTFail,-deltaTFail)        #curTime is in seconds
#     # timeVec = np.arange(tfail*1.0,0,-deltaTFail)        #curTime is in seconds
#
#     jobs = [job_server.submit(MonteCarlo_Distributed_Reentry_and_Record, \
#                               args=(curMission, desiredCatalog, coeffIX, numPiecesPerSample, \
#                                                  profiles, lowerBreakLimit, upperBreakLimit, windIX), \
#                               depfuncs=(MonteCarlo_Distributed_Reentry,), \
#                               modules=('numpy as np','debrisReader as DR', 'orbitTools', 'debrisPropagation as dp'), \
#                               callback=finished001) for windIX in range(numWindSamples)]
#
#     # jobs = [job_server.submit(MonteCarlo_at_tfail_and_record, \
#     #                           args=(curMission, coeffIX, curTime, numPiecesPerSample, profiles), \
#     #                           depfuncs=(MonteCarlo_at_tfail,), \
#     #                           modules=('numpy as np','debrisReader as DR', 'orbitTools', 'debrisPropagation as dp'), \
#     #                           callback=finished001) for curTime in timeVec]
#
#     job_server.wait()
#     job_server.print_stats()
#     for job in jobs:
#         print job()
#
#     # return cur_mpc
#
#
#
#
#
# ### Wrapper function used for efficient parallelization of MonteCarlo_until_tfail (Called by MonteCarlo_until_tfail)
# def MonteCarlo_Distributed_Reentry_and_Record(curMission, desiredCatalog, coeffIX, numPiecesPerSample,
#                                                  profiles, lowerBreakLimit, upperBreakLimit, windIX):
#
#     # mpc = MonteCarlo_at_tfail(curMission, coeffIX, curTime, numPiecesPerSample, profiles)
#
#     cur_mpc = MonteCarlo_Distributed_Reentry(curMission, desiredCatalog, coeffIX, numPiecesPerSample,
#                                                  profiles, lowerBreakLimit, upperBreakLimit, windIX)
#
#     debrisPickleFolder = curMission['debrisPickleFolder']
#
#     # print 'COMMENTED OUT WRITING TO FILE FOR DEBUGGING PURPOSES'
#     # Make sure that the output directory exists
#     folderPath = os.path.abspath(debrisPickleFolder)
#     if not os.path.exists(folderPath):
#         os.makedirs(folderPath)
#
#     output = open(folderPath + '/mpc_' + str(windIX) + '.pkl', 'wb')
#     pickle.dump(cur_mpc,output,2)
#     output.close()
#
#     return windIX
#
#     # return curTime
#
#
# ### CURRENTLY IN USE (Called by Columbia.py -->
# def MonteCarlo_Distributed_Reentry(curMission, desiredCatalog, coeffIX, numPiecesPerSample, profiles, lowerBreakLimit, upperBreakLimit, windIX):
#
#     debrisCatalogFilePATH = curMission['debrisCatPath']
#     dt = curMission['deltaT']
#
#     ndtinterval = 5*3600/dt     # is the upper bound for how long you think the debris propagation will run (5hrs)
#     tLaunchDesired = 0.         # Delay in seconds of launch from nominal time
#
#     # Unpack the profiles
#     atmStorage = profiles['atmStorage']
#     thetagStorage = profiles['thetagStorage']
#     stateVecStorage = profiles['stateVecStorage']
#     tfailStorage = profiles['tfailStorage']
#
#     # Find sizes of profiles
#     # numWindSamples = len(atmStorage)
#     numTrajSamples = len(stateVecStorage[0])
#
#     omegaE = curMission['omegaE']
#     planetmodel = curMission['planetModel']
#
#     # These two arrays will be appended to and returned as the solution of this function
#     debrisNumTimeSteps = []
#     debrisID = []
#     debrisMass = []
#     debrisArea = []
#     testStorage = []
#
#     altitudeList    = atmStorage[windIX][0]
#     densityList     = atmStorage[windIX][1]
#     uList           = atmStorage[windIX][2]
#     vList           = atmStorage[windIX][3]
#     wList           = atmStorage[windIX][4]
#
#     for trajIX in range(numTrajSamples):
#
#         # generating random debris pieces
#         massList,arefList,velList,ncdList,minfcdList,CDList,nclList,minfclList,CLList,arefMean,numberOfPiecesMeanList,blast=DR.generateDebris(desiredCatalog,numPiecesPerSample,debrisCatalogFilePATH)
#
#
#         for debrisIndex in coeffIX:
#             ncdval = ncdList[debrisIndex]
#             minfcdval = minfcdList[debrisIndex]
#             nclval = nclList[debrisIndex]
#             minfclval = minfclList[debrisIndex]
#             CDval = CDList[debrisIndex]
#             CLval = CLList[debrisIndex]
#
#             # print '[wind][traj][beta] = ' + '[{0}][{1}][{2}]'.format(windIX, trajIX, debrisIndex)
#
#             for pieceIndex in range(numPiecesPerSample):
#
#                 # sampling time during a two minute window
#                 tbreak = np.random.uniform(lowerBreakLimit,upperBreakLimit)
#
#                 # This gets moved into loop for each piece within a generated debris catalog
#                 nearestIX = (np.abs(tfailStorage[0][0]-tbreak)).argmin() #Assuming first windprofile first trajectory is representative
#
#                 thetag = thetagStorage[windIX][trajIX][nearestIX]
#                 x = stateVecStorage[windIX][trajIX][0][nearestIX]
#                 y = stateVecStorage[windIX][trajIX][1][nearestIX]
#                 z = stateVecStorage[windIX][trajIX][2][nearestIX]
#                 Vx = stateVecStorage[windIX][trajIX][3][nearestIX]
#                 Vy = stateVecStorage[windIX][trajIX][4][nearestIX]
#                 Vz = stateVecStorage[windIX][trajIX][5][nearestIX]
#
#                 stateVec = [x,y,z,Vx,Vy,Vz]
#
#                 Vrot = np.cross([0,0,omegaE],[x,y,z])
#                 Vinf = np.array([Vx,Vy,Vz]) - Vrot
#                 Vmag = (Vinf[0]**2 + Vinf[1]**2 + Vinf[2]**2)**.5
#
#                 massval = massList[debrisIndex][pieceIndex]
#                 arefval = arefList[debrisIndex][pieceIndex]
#                 velImpval = velList[debrisIndex][pieceIndex]
#                 #velImpval = 100.  # addition...just to check behavior
#                 theta1 = np.random.uniform(0.0,2*np.pi)
#                 theta2 = np.random.uniform(0.0,2*np.pi) #it is not required this to be 0 to 2pi, domain is covered with 0 to pi since prev angle helps cover the entire sphere
#                 Vdeb = np.array([[0.0*velImpval,0,0]]).T # velocity impulse calculation
#                 Rz = orbitTools.cRnMatrix_Zaxis(theta1)
#                 Rx = orbitTools.cRnMatrix_Xaxis(theta2)
#                 Vdeb = np.dot(Rz,Vdeb)
#                 Vdeb = np.dot(Rx,Vdeb)
#                 Vdeb = Vdeb[:,0] # making it a 1 by 3 vector...impulse velocity in cartesin coordinates
#                 #print arefval
#                 # debris propagation calculation
#
#                 # print 'pieceIndex = {0}, Vdeb/Vmag = {1}/{2}'.format(pieceIndex, (Vdeb[0]**2 + Vdeb[1]**2 + Vdeb[2]**2)**.5, Vmag)
#
#                 #            print 'pieceIX = ' + str(pieceIndex) + '   mass = ' + str(massval)
#
#                 debrisResults, numFinalSteps = dp.debrispropagation(initialstate = stateVec,
#                                                                     debrisvel = Vdeb,
#                                                                     mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=1,
#                                                                     minfcl=minfclval,cl=CLval,loverd = 0.04, atmosoption=2,altitudelist=altitudeList,
#                                                                     densitylist=densityList,ulist=uList,vlist=vList,
#                                                                     wlist=wList,geoptions=0,filename='none',planetmodel=planetmodel,dtinterval = dt,ndtinterval = ndtinterval,thetag0=thetag,
#                                                                     ncd=ncdval,
#                                                                     ncl=nclval,
#                                                                     nlist=len(altitudeList))
#                 #                pdb.set_trace()
#
#     #                    print 'QUITTING AFTER ONE ITERATION AND RETURNING TEST VALUE FOR DEBUGGING PURPOSES!!!!'
#     #                    return debrisResults, numFinalSteps
#
#                 # Explanation of these manipulations:
#                 # debrisResults is returned as a fixed-length one-dimensional array.  Fortran makes it way larger than it
#                 #   would ever need to be because it doesn't have the ability to do dynamic memory allocation and it doesn't
#                 #   know at the beginning of the propagation how long it will take the debris to hit the ground.  Thus, we need
#                 #   to trim all the unused stuff off of the end
#                 #   * debrisResults[0:(numFinalSteps)*sizeResults]
#                 # Then we want to unflatten the array into a multidimensional array
#                 #   * .reshape(numFinalSteps,sizeResults)
#                 # But there might be return values that we don't care about, so select only the columns that we want to store.
#                 #   In my case, that's the first six values (position and velocity)
#                 #   * [:,:6]
#                 # Later when we store the values in debrisStorage, we'll flatten them out and make this is a massive 1D array
#                 #   * .flatten()
#                 sizeResults = 7
#                 sizeKeeping = 6
#                 debrisResults = debrisResults[0:(numFinalSteps)*sizeResults].reshape(numFinalSteps,sizeResults)[:,:sizeKeeping]
#
#                 altitudeFinal = debrisResults[-1,2]
#                 if altitudeFinal>0:
#                     print 'Warning...debris is not on the ground. Final altitude is ',altitudeFinal
#
#                 # Saves the stuff I care about
#                 # For some reason, those extra parentheses are necessary to avoid errors when concatenating the first (empty) array
#                 # debrisStorage = np.concatenate(( debrisStorage,  debrisResults.flatten() ))
#                 testStorage.append(debrisResults.flatten() )
#                 debrisNumTimeSteps.append(numFinalSteps)
#                 debrisID.append(debrisIndex)
#                 debrisMass.append(massval)
#                 debrisArea.append(arefval)
#
#     # Convert the time steps array to be numpy int32 array
#     debrisNumTimeSteps = np.array(debrisNumTimeSteps,dtype=np.int32)
#     debrisID = np.array(debrisID,dtype=np.int32)
#     debrisMass = np.array(debrisMass,dtype=np.double)
#     debrisArea = np.array(debrisArea,dtype=np.double)
#     numberOfPiecesMeanList = np.array(numberOfPiecesMeanList,dtype=np.int32)
#     arefMean = np.array(arefMean,dtype=np.double)
#
#     # This flattens the python list testStorage and saves it as the proper numpy array
#     debrisStorage = np.array([item for sublist in testStorage for item in sublist],dtype=np.float64)
#
#
#     ###### Now to read these points into a PointCloud, for input into a footprint and export to GE
#     # UTC = curMission['initialUTC'] + (tLaunchDesired + tfail)/(24*3600.)
#     UTC = curMission['initialUTC']
#     # zBinHeightKm = 5.0
#
#     numRuns = len(debrisNumTimeSteps)
#     maxTime = debrisNumTimeSteps.max()*dt
#
#     # Pack things up in a dictionary (for pickling)
#     mpc = dict(flatPointArray = debrisStorage, numPieces = numRuns, numTimeSteps = debrisNumTimeSteps, maxTime = maxTime, deltaTsec = curMission['deltaT'],
#                UTC = UTC, all_points_delta_t = curMission['all_points_delta_t'],
#                launchLat = curMission['launchLat'], launchLon = curMission['launchLon'], launchAzimuth = curMission['launchAzimuth'],
#                debrisID = debrisID, arefMeanList = arefMean, numberOfPiecesMeanList = numberOfPiecesMeanList,
#                debrisMass = debrisMass, debrisArea = debrisArea, sizeFlatPointArray = sizeKeeping)
#
#
#
#     return mpc






















# def InvestigateDebrisCatalog(debrisCatalogFile):
#
#     # getting debris catalog information, getting desired debris catalog group
#     velcat,timecat,catalogList = DR.readDebris(debrisCatalogFile)
#     desiredCatalog = DR.desiredList(velcat,timecat,catalogList,Vmag,tfail)
#
#     # generating random debris pieces
#     massList,arefList,velList,ncdList,minfcdList,CDList,nclList,minfclList,CLList,arefMean,numberOfPiecesMeanList,blast=DR.generateDebris(desiredCatalog,numPiecesPerSample,debrisCatalogFilePATH)






def InitializeMission(mission1):
    # Unpack the dictionary
    propagationParamFile    = mission1['propagationParamFile']
    precomputedParamFile    = mission1['precomputedParamFile']
    pathToMissionFiles      = mission1['pathToMissionFiles']

    # Load the appropriate files
    if len(propagationParamFile) > 0:
        # This is a trajectory that we need to propagate.  Load up the parameters
        from ReadMission import readInput
        newMission = readInput(propagationParamFile, pathToMissionFiles)
    elif len(precomputedParamFile) > 0:
        # This is a precomputed trajectory, load up the state vectors as best you can
        newMission = LoadPrecomputed(mission1)
        # newMission['atmoOption'] = 1
    else:
        print 'Something went terribly wrong.  No mission loaded'
        newMission = dict()

    # Place the inputs from the original mission into the dict as well
    newMission['propagationParamFile']  = propagationParamFile
    newMission['precomputedParamFile']  = precomputedParamFile
    newMission['pathToMissionFiles']    = pathToMissionFiles
    newMission['pathToMissionFiles']    = pathToMissionFiles
    newMission['omegaE']                = mission1['omegaE']
    newMission['planetModel']           = mission1['planetModel']
    return newMission



def LoadPrecomputed(mission):
    # Load the latlonalt profile
    # FlightProfileFile = pathToMissionFiles + 'HTHLwCarrier.txt'
    omegaE = mission['omegaE']
    planetModel = mission['planetModel']

    maxVRel = 0.;
    inputFile = open(mission['pathToMissionFiles'] + mission['precomputedParamFile'], 'r')

    stateVecProfile = []
    from FriscoLegacy import orbitTools as ot

    isLATLON    = False
    isSTATEVEC  = False

    isMETERS    = False
    ft2m = 0.3048

    firstLine = inputFile.readline()
    key = firstLine.split()
    if (key[0] == 'FORMAT'):
        if (key[1] == 'LATLON'):
            isLATLON = True
            if (len(key) == 3) and (key[2].upper() == 'METERS'):
                isMETERS = True
                ft2m = 1.   # kind of a hack
            elif (len(key) == 3) and (key[2].upper() == 'FEET'):
                isMETERS = False
                # ft2m = 1.   # kind of a hack
            else:
                print 'ERROR! You did not specify the units in your precomputed file.  Must be METERS or FEET'
                raise RuntimeError

        elif (key[1] == 'STATEVEC'):
            isSTATEVEC = True
        else:
            print key[1] + ' is not a recognized format.  Must be either LATLON or STATEVEC.'
            raise RuntimeError
    else:
        print 'first line keys = ' + str(key)
        print 'First line of file must read: FORMAT LATLON or FORMAT STATEVEC\n\n'
        raise RuntimeError


    # print 'inputFile[0] = ' + inputFile.readline()

    if isLATLON:

        tempStorage = []    # Used for holding on to info from first time step
        count = -1

        for line in inputFile:
            key = line.split()

            #if line is not empty and not a comment
            if (len(key)>0) and (key[0][0] != '#'):
                count = count + 1   # Count starts at -1, so first time through this loop count will be 0

                # Altitude units are feet???
                ft2km = ft2m * (1e-3)

                # Convert the strings into floats
                curTime = float(key[0])
                curLat = float(key[2])
                curLon = float(key[1])
                curAltft = float(key[3])
                # timelatlonaltProfile.append((float(key[0]), float(key[2]), float(key[1]), curAltft*ft2km))

                # Basically just use the start of the mission as the reference time and don't worry about greenwich
                thetag = omegaE * curTime

                # Attempt conversion right here
                [x_eci,y_eci,z_eci] = ot.latlonalt2ECI(curLat, curLon, curAltft*ft2m, thetag, planetModel)


                # Simple backward difference
                if count > 1:
                    # This is not one of the first two timesteps, so nothing special needs to happen
                    delta_t = curTime - stateVecProfile[-1][0]
                    Vx = (x_eci - stateVecProfile[-1][2])/delta_t
                    Vy = (y_eci - stateVecProfile[-1][3])/delta_t
                    Vz = (z_eci - stateVecProfile[-1][4])/delta_t

                    Vrotfinal = np.array([-y_eci*omegaE, x_eci*omegaE, 0.])
                    VrelMag = np.linalg.norm(np.array([Vx, Vy, Vz])-Vrotfinal)

                    stateVecProfile.append((curTime, thetag, x_eci, y_eci, z_eci, Vx, Vy, Vz, VrelMag))

                    # maxVRel = np.max( (maxVRel,VrelMag) )
                    # print 'VrelMag = ' + str(VrelMag*(1e-3))
                    # print str(curTime) + ' VrelMag = ' + str(VrelMag)

                elif count == 0:
                    # This is the first time through the loop and we don't have enough info yet to do a backwards difference
                    # So just save the xyz portion of the state vector for later
                    tempStorage = [curTime, thetag, x_eci, y_eci, z_eci]

                elif count == 1:
                    # This is the second time through, so you need to handle the rest of the first timestep
                    # Unpack the previously stored values
                    prevCurTime = tempStorage[0]
                    prevThetag  = tempStorage[1]
                    prevX_eci   = tempStorage[2]
                    prevY_eci   = tempStorage[3]
                    prevZ_eci   = tempStorage[4]

                    # Now, use that to figure out the velocities (Same for both time steps)
                    delta_t = curTime - prevCurTime
                    Vx = (x_eci - prevX_eci)/delta_t
                    Vy = (y_eci - prevY_eci)/delta_t
                    Vz = (z_eci - prevZ_eci)/delta_t

                    Vrotfinal = np.array([-y_eci*omegaE, x_eci*omegaE, 0.])
                    VrelMag = np.linalg.norm(np.array([Vx, Vy, Vz])-Vrotfinal)

                    # Push them both, in order, into the state vector storage
                    stateVecProfile.append((prevCurTime, prevThetag, prevX_eci, prevY_eci, prevZ_eci, Vx, Vy, Vz, VrelMag))
                    stateVecProfile.append((curTime, thetag, x_eci, y_eci, z_eci, Vx, Vy, Vz, VrelMag))

                    # Destory the temp
                    tempStorage = []


    elif isSTATEVEC:
        for line in inputFile:
            key = line.split()

            #if line is not empty and not a comment
            if (len(key)>0) and (key[0][0] != '#'):
                curTime     = float(key[0])
                thetag      = float(key[1])
                x_eci       = float(key[2])
                y_eci       = float(key[3])
                z_eci       = float(key[4])
                Vx          = float(key[5])
                Vy          = float(key[6])
                Vz          = float(key[7])

                VrelMag     = 1     # This doesn't really get used, so whatever
                stateVecProfile.append((curTime, thetag, x_eci, y_eci, z_eci, Vx, Vy, Vz, VrelMag))



    # Convert to numpy array
    stateVecProfile = np.array(stateVecProfile)

    # # Needed (?)
    # # These should probably get set in a param file eventually
    # atmoOption = 1.
    # atmoFile = 'special.txt'
    # mission1 = dict(atmoOption = atmoOption, atmoFile = atmoFile, stateVecProfile = stateVecProfile)

    mission1 = dict(stateVecProfile = stateVecProfile)

    return mission1





# propagationParamFile = []
# precomputedParamFile = 'HTHLwCarrier.txt'






















# # NOT CURRENTLY BEING USED.  However, I think I will want this in the future.
# def GenerateLRHC(curMission, numRuns):
#     print 'In GenerateLRCH'
#
#     # getting atmospheric profile
#     altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints= AP.readGramAtmos(curMission['atmosphere'])
#     # making sure atmospheric profile is decreasing in altitude
#     altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist= AP.maxAltitudeAtmosFix(altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints)
#
#     from scipy.io import loadmat
#     import orbitTools as ot
#
#     missionList = loadmat(curMission['trajMatFile'])['missionList'] # loading trajectories  from SPOT
#     trajStorage = np.array([])
#     trajNumTimeSteps = []
#     maxTime = 0;
#
#     for tx in range(numRuns):
#         print 'tx = ' + str(tx)
#
#         tLaunchDesired = 0.0 # launch time in seconds from beginning of launch window
#         ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#         Toffset = np.random.normal(1,.1) # value to multiply the thrust...
#         densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList) # generating new atmos profile
#         atmProfile1 = [altitudeList,densityList,uList,vList,wList]
#
#         myTrajList = traj.generateRandomTrajFAST(missionList,tLaunchDesired,curMission['deltaT'],atmProfile1,Toffset,ThrustOffsetAngDeg)
#
# #        myTrajList = traj.generateRandomTraj(missionList,tLaunchDesired,curMission['deltaT'],atmProfile1,Toffset,ThrustOffsetAngDeg)
#         [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent] = myTrajList[0]    #index 0=first stage, 1=coast period, 2=second stage
#
#         rvec = np.array((xc, yc, zc))
#         planetModel = 0     #spherical earth
#
#
#         latlonalt = np.zeros((3,xc.size))
#
#         for ix in range(0,xc.size):
#             latlonalt[:,ix] = ot.ECI2latlonalt(rvec[:,ix],thetagVector[ix],planetModel)
#
#         trajStorage = np.concatenate((trajStorage,latlonalt.flatten(1)))    #later just load these in properly, but for now flatten them into a 1d array like debris
#         #    trajNumTimeSteps = np.concatenate((trajNumTimeSteps, tc.size))
#         trajNumTimeSteps.append(tc.size)
#         if (tc.max() > maxTime):
#             maxTime = tc.max()
#
#     trajNumTimeSteps = np.array(trajNumTimeSteps,dtype=np.int32)    #convert into int32 array
#
#     ###### Now to read these points into a PointCloud, for input into a footprint and export to GE
#     UTC = 156.848611111
#     zBinHeightKm = 5.0
#
#     launchLat = 28.445455
#     launchLon = -80.564865
#     launchAzimuth = 0.976192685683
#
#     all_points_delta_t = 1.0    #minutes
#
#     #(&flatPointArray[0], numPieces, &numTimeSteps[0], maxTime, deltaTsec)
#     myPointCloud = ceb.PyPointCloud(trajStorage, numRuns, trajNumTimeSteps, maxTime, curMission['deltaT'],
#                                     UTC, all_points_delta_t,
#                                     launchLat, launchLon, launchAzimuth)
#
#     return myPointCloud













# # SEEMS UNUSED!!!
# def generateRandomTraj(missionList,tLaunchDesired,dt,atmProfile=[],TOffset = 1.0,ThrustOffsetAngDeg=[0,0]):
#     # this function generates entire trajectories
#     import numpy as np
#     #calculating ECI coordinates for launch
#     import orbitTools # library for coordinate transformation
#     import orbitProp as op # trajectory propagation routine
#     import data2GE
#     import copy
#
#     if len(atmProfile)>0:
#         atmosoption =2
#         [altitudelist,densitylist,ulist,vlist,wlist] = atmProfile
#     else:
#         atmosoption = 0
#         densitylist = [0]
#         ulist = [0]
#         vlist = [0]
#         wlist = [0]
#         altitudelist = [0]
#
#     mission0 = missionList[0][0][0]# getting mission parameters for current mission
#     omegaE = mission0['planet'][0][0]['omega'][0][0][0] # planet's angular velocity
#     retList = []
#
#     latitude = mission0['launchSite'][0][0]['lat'][0][0][0]
#     longitude = mission0['launchSite'][0][0]['long'][0][0][0]
#     height = mission0['launchSite'][0][0]['h'][0][0][0]
#
#     planetmodel = 0# 0 for Circ ...1 for ellipt
#     trajList,thetag0 = interpolateTrajectories(missionList,tLaunchDesired,dt)
#
#     r0 = orbitTools.latlonalt2ECI(latitude,longitude,height,thetag0,planetmodel)
#
#     nEvents = mission0['vehicle'][0][0]['stages'][0][0][0]
#     cloption = 1
#     cl = 0.
#     minfcl = [1]
#     loverd = 0.0
#
#     geoptions = 0
#
#     dtinterval = dt
#
#     ndtinterval = 20000
#     #    ndtinterval = 60
#     fileList = []
#     currentTime = 0.0
#     retList= []
#     for indexEvent in range(nEvents):
#         if indexEvent ==0 : #first stage...setting initial conditions for rest of cases
#
#             [vS0,vE0,vZ0] = [0,0,0] # initial velocity in SEZ frame
#             inertial = 0
#             thetag = copy.deepcopy(thetag0)
#             # inertial = 1. input velocities are inertial velocities in SEZ frame
#             # inertial =0. input velocities are local velocities (earth rotation not accounted)
#             #Veci0 = orbitTools.SEZ2ECI(latitude,longitude,r0,vS0,vE0,vZ0,inertial,thetag)
#             Veci0 = orbitTools.SEZ2ECI(latitude,longitude,height,vS0,vE0,vZ0,inertial,thetag,planetmodel)
#             initialstate = np.array([r0[0],r0[1],r0[2],Veci0[0],Veci0[1],Veci0[2]])
#
#         sref =  mission0['vehicle'][0][0]['aero'][0][0]['Aref'][0][indexEvent][0]
#         cd = mission0['vehicle'][0][0]['aero'][0][0]['CD'][0][indexEvent][0][:,0]
#         minfcd = mission0['vehicle'][0][0]['aero'][0][0]['MinfCD'][0][indexEvent][0][:,0]
#         mass0 = mission0['vehicle'][0][0]['mass'][0][0]['m0'][0][indexEvent][0]
#         massf = mission0['vehicle'][0][0]['mass'][0][0]['mf'][0][indexEvent][0]
#         isp =  mission0['vehicle'][0][0]['propulsion'][0][0]['ISP'][0][indexEvent][0]
#         # Tmax including constant offset for all stages
#         Tmax = TOffset*mission0['vehicle'][0][0]['propulsion'][0][0]['Ftotal'][0][indexEvent][0]
#         etaMin = mission0['vehicle'][0][0]['propulsion'][0][0]['minThrottle'][0][indexEvent][0]
#         filename = 'pTraj'+str(indexEvent)+'.txt'
#         fileList.append(filename)
#         if etaMin>=0:
#
#             propOption = 1
#             currTraj = trajList[indexEvent]
#             [indexEvent,tfit,uxFit,uyFit,uzFit,etaFit] = currTraj
#             # ensuring that unit vector is actually a unit vector (due to interpolation)
#
#             uMag = (uxFit**2 + uyFit**2 + uzFit**2)**.5
#             uxetaFit = np.array([etaFit*uxFit/uMag]).T
#             uyetaFit = np.array([etaFit*uyFit/uMag]).T
#             uzetaFit = np.array([etaFit*uzFit/uMag]).T
#             uetaMat = np.concatenate((uxetaFit,uyetaFit,uzetaFit),1)
#             Tmat = Tmax*uetaMat
#             timelist = tfit
#             propCond = massf
#             ntime = len(timelist)
#
#             Rearth = mission0['planet'][0][0]['R'][0][0][0]
#             omega = mission0['planet'][0][0]['omega'][0][0][0]
#
#
#             finalconditions,finalderivs = op.propagate(initialstate, mass0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Tmat,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
#             newIndex = finalconditions[:,10]>=0.0
#             newfinalConditions = finalconditions[newIndex,:]
#             # updating parameters for propagation in next stage
#             xc = newfinalConditions[:,0]
#             yc = newfinalConditions[:,1]
#             zc = newfinalConditions[:,2]
#             Vxc = newfinalConditions[:,3]
#             Vyc = newfinalConditions[:,4]
#             Vzc = newfinalConditions[:,5]
#             mc = newfinalConditions[:,9]
#             tc = newfinalConditions[:,10] + currentTime
#             initialstate = np.array([xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]])
#             mass0 = mc[-1]
#
#             currentTime = tc[-1] # updating time
#             thetagVector = omegaE*tc + thetag0
#             thetag = omegaE*currentTime + thetag0 # updating thetag for next event
#
#         elif etaMin<0:
#
#             propOption = 2
#             currTraj = trajList[indexEvent]
#             [indexEvent,tffit] = currTraj
#             dtlocal = tffit/5.
#             Tmax = 0.0
#             uetaMat = np.array([[1,1,1]])
#             Tmat = Tmax*uetaMat
#             timelist = [tffit]
#             propCond = tffit
#             ntime = len(timelist)
#
#
#
#             finalconditions,finalderivs = op.propagate(initialstate, mass0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Tmat,timelist,isp,geoptions,filename,planetmodel,dtlocal,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
#             newIndex = finalconditions[:,10]>=0.0
#             newfinalConditions = finalconditions[newIndex,:]
#             # updating parameters for propagation in next stage
#             xc = newfinalConditions[:,0]
#             yc = newfinalConditions[:,1]
#             zc = newfinalConditions[:,2]
#             Vxc = newfinalConditions[:,3]
#             Vyc = newfinalConditions[:,4]
#             Vzc = newfinalConditions[:,5]
#             mc = newfinalConditions[:,9]
#             tc = newfinalConditions[:,10] + currentTime
#             initialstate = np.array([xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]])
#             mass0 = mc[-1]
#             currentTime = tc[-1] # updating time
#             thetagVector = omegaE*tc + thetag0
#             thetag = omegaE*currentTime + thetag0# updating thetag for next event
#
#         retList.append([xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent])
#     return retList
# #data2GE.convertMultiple(fileList)







# Initial velocity of zero is hardcoded; should read in initial velocity
def getStateVector(curMission, atmProfile, Toffset, ThrustOffsetAngDeg):
    # import orbitTools
    from scipy.interpolate import UnivariateSpline
    # from scipy import interpolate
    from FriscoLegacy import orbitProp as op
    import sys
    sys.path.append('../PythonScripts')
    # import AtmosProfile as AP
    #import matplotlib.pyplot as plt
    planetModel = 0
    
    #outputs = [initialLocation,massList,ArefList,MinfCDlist,CDList,atmoOption,atmoFile,timelist,thrustList,VehicleIsp,thetag0,dateList,dtval]
    
    omegaE = 2.*np.pi/(86164.0906)
    
    #print curMission.keys()
    # Unpack the trajectory information
    [lon0,lat0,h0]  = curMission['initialLocation']
    massList        = curMission['massList']
    SrefList        = curMission['ArefList']
    MinfCDlist      = curMission['MinfCDlist']
    CDList          = curMission['CDList']
    atmosoption     = curMission['atmoOption']
    atmoFile        = curMission['atmoFile']
    timelist        = curMission['timelist']
    thrustList      = curMission['thrustList']
    ispList         = curMission['VehicleIsp']
    thetag          = curMission['thetag0']
    # Note that curMission['dateList'] appears to be unused
    # dtinterval      = curMission['dtval']
    dtinterval      = curMission['deltaT']
    ndtinterval     = curMission['debrisTimeLimitSec']/dtinterval


    cloption = 1
    if curMission['useLoverD'] == True:
        cloption = 0
    loverd = curMission['loverd']


    [x0,y0,z0] = orbitTools.latlonalt2ECI(lat0,lon0,h0,thetag,planetModel)
    #SEZ2ECI(lat,lon,hs,Vsouth,Veast,Vzenith,inertial,thetag,planetModel)
    V0 = orbitTools.SEZ2ECI(lat0,lon0,h0,0,0,0,0,thetag,planetModel)
    Vx0 = V0[0]
    Vy0 = V0[1]
    Vz0 = V0[2]
    m0Vec = massList[0]
    mfVec = massList[1]

    retTrajs = []       # This will store the trajectories that get returned
    #atmosoption = 2
    if atmosoption>0:
        print 'using given atmo data to propagate a trajectory'
        #atmStorage.append([altitudeList,densityList,uList,vList,wList])
#        atmosoption =2
        [altitudeList,densitylist,ulist,vlist,wlist] = atmProfile
#        altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints= AP.readGramAtmos(atmoFile)
    else:
        altitudeList = [0]
        densityMeanList = [0]
        uMeanList = [0]
        vMeanList = [0]
        wMeanList = [0]
    
#    densitylist = densityMeanList
#    ulist = uMeanList
#    vlist = vMeanList
#    wlist = wMeanList
    
    # cloption = 0
    cl = [0]
    # loverd = 0
    minfcl = [1]
    geoptions = 1
    fileGE = curMission['GeneratedFilesFolder'] + 'trajProp'
#    thrustOffsetAngDeg =[0,0]
    mass_time_alt_opt = 1
#    dtinterval = 0.01
    
    currentTime = 0;
    
    stateVectors = []
    for index in range(0,len(m0Vec)):
        mass    = m0Vec[index]
        mf      = mfVec[index]
        thrust  = thrustList[index]
        timeVec = timelist[index]
        cd0     = CDList[index]
        minfcd0 = MinfCDlist[index]
        sref    = SrefList[index]
        isp     = ispList[index]
        
        initialstate = [x0,y0,z0,Vx0,Vy0,Vz0]   #These get updated at end of loop

#        Tx = thrust[0]
#        Ty = thrust[1]
#        Tz = thrust[2]
        
#        # Assume that halfway through the stage, we're at full thrust
#        Tmax = np.linalg.norm( np.transpose(thrust)[np.floor(len(Tx)/2),:])
#        
#        # But the old Tmax was actually the maximum OFFSET, so...
#        Tmax = Tmax * Toffset
        
        # May not even need Tmax...i think this is the way to handle it
        Tx = thrust[0] * Toffset
        Ty = thrust[1] * Toffset
        Tz = thrust[2] * Toffset
        
        
        '''
            fTx = UnivariateSpline(timeVec,Tx)
            fTy = UnivariateSpline(timeVec,Ty)
            fTz = UnivariateSpline(timeVec,Tz)
            
            timeVec = np.linspace(timeVec[0],timeVec[-1],1000)
            Tx = fTx(timeVec)
            Ty = fTy(timeVec)
            Tz = fTz(timeVec)
            '''
        
        fcd0 = UnivariateSpline(minfcd0,cd0)
        minfcd = np.linspace(minfcd0[0],minfcd0[-1],100)
        cd = fcd0(minfcd)
        Tlist = np.array([Tx,Ty,Tz]).T#np.concatenate([[Tx],[Ty],[Tz]],axis=1)
        
        
        # ndtinterval = int(np.ceil(2*timeVec[-1]/dtinterval))
        ntime = len(timeVec)
        #finalconditions,derivVals = op.propagate(initialstate,mass,mf,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudeList,densitylist,ulist,vlist,wlist,Tlist,timeVec,isp,geoptions,fileGE,planetModel,dtinterval,ndtinterval,thetag,1,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudeList))
        
        #finalconditions,finalderivs = propagate(initialstate,mass,mass_time_alt_final,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,tlist,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,mass_time_alt_opt,thrustoffangledeg,ncd=len(minfcd),ncl=len(minfcl),ntime=shape(tlist,0),nlist=len(altitudelist))
        

        
        if index == len(m0Vec)-1:
            mass_time_alt_opt = 2
            fcond = timeVec[-1]
        else:
            fcond = mf
#        # New
#        finalconditions,finalderivs = op.propagate(initialstate,mass,fcond,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudeList,densitylist,ulist,vlist,wlist,Tlist,timeVec,isp,geoptions,fileGE+str(1+index),planetModel,dtinterval,ndtinterval,thetag,mass_time_alt_opt)

        # Incorporating the Old
        finalconditions,finalderivs = op.propagate(initialstate,mass,fcond,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudeList,densitylist,ulist,vlist,wlist,Tlist,timeVec,isp,geoptions,fileGE+str(1+index),planetModel,dtinterval,ndtinterval,thetag,mass_time_alt_opt,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudeList))
        
        newIndex = finalconditions[:,10]>=0.0
        newfinalConditions = finalconditions[newIndex,:]
        stateVectors.append(newfinalConditions)
        
        x0 = newfinalConditions[-1,0]
        y0 = newfinalConditions[-1,1]
        z0 = newfinalConditions[-1,2]
        Vx0 = newfinalConditions[-1,3]
        Vy0 = newfinalConditions[-1,4]
        Vz0 = newfinalConditions[-1,5]
        massCheck = newfinalConditions[-1,9]
        timeCheck = newfinalConditions[-1,7]
        
        #thetag = thetag + omegaE*newfinalConditions[-1,10]
        
        # Format things for output (Just get it to work first)
        xc = newfinalConditions[:,0]
        yc = newfinalConditions[:,1]
        zc = newfinalConditions[:,2]
        Vxc = newfinalConditions[:,3]
        Vyc = newfinalConditions[:,4]
        Vzc = newfinalConditions[:,5]
        mc = newfinalConditions[:,9]
        tc = newfinalConditions[:,10] + currentTime
        
        thetagVector = curMission['thetag0'] + omegaE*tc
        thetag = thetagVector[-1]
        
        currentTime = tc[-1] # updating time
        
        
#        mass0 = mc[-1]
#        thetagVector = omegaE*tc + thetag0
#        thetag = omegaE*currentTime + thetag0 # updating thetag for next event
        
        retTrajs.append([xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,index])
    return retTrajs












# ========================= GRAVEYARD =========================== #


#def ensure_dir(f):
#    d = os.path.abspath(f)
#    if not os.path.exists(d):
#        os.makedirs(d)







##(atmoFile, trajectoryFiles, debrisCatalogFile, debrisCatalogFilePATH, dt, tLaunchDesired, tfail, nSamples, nRuns):
#def MonteCarloSingleBallisticCoeff(curMission, coeffIX, tfail, numTrajSamples, numWindSamples, numPiecesPerSample):
#    print 'In MonteCarloSingleBallisticCoeff, tLaunchDesired = ZERO'
#
#    #    dt = 2.0 # stepping for trajectory propagation
#
#    # Unpack some of the mission variables
#    tLaunchDesired = 0.
#
#    debrisCatalogFile = curMission['debrisCatPath'] + curMission['debrisCatFile']
#    debrisCatalogFilePATH = curMission['debrisCatPath']
#    atmoFile = curMission['atmosphere']
#    dt = curMission['deltaT']
#
#    tfailIX = tfail/dt;
#    if (np.mod(tfail,dt) != 0.0):
#        print 'YOU PICKED A BAD FAILURE TIMESTEP'
#        return 0
#
#
#
#
#    ndtinterval = 5*3600/dt    #is the upper bound for how long you think the debris propagation will run (5hrs)
#
#    # getting atmospheric profile
#    altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints= AP.readGramAtmos(atmoFile)
#    # making sure atmospheric profile is decreasing in altitude
#    altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist= AP.maxAltitudeAtmosFix(altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints)
#
#
#
#    # getting trajectories
#    #trajVals = traj.generateTraj(trajectoryFiles,atmProfile,5)
#    #stateVectors,time,thetag = traj.pickTraj(trajVals,0,15,20)
#    #retList = traj.getALLTraj(trajectoryFiles,.1)
#
#    #getting tf
#    missionList = loadmat(curMission['trajMatFile'])['missionList'] # loading trajectories  from SPOT
#
#
##    import time
#
#
#    #######debris catalog selection and generation#################################
#    # Need a trajectory first, which requires a wind profile, in order to get a Vmag.
#    # Make temporary ones of both and have them fall out of scope
#    Vmag = -5.
#
#    Toffset_std = 0.02;
#
#    if (1 == 1):
#        densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
#        atmProfile1 = [altitudeList,densityList,uList,vList,wList]
#        ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#        Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
#        stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile1,Toffset,ThrustOffsetAngDeg)
#
#    print 'Vmag = ' + str(Vmag)
#
#    # getting debris catalog information
#    velcat,timecat,catalogList = DR.readDebris(debrisCatalogFile)
#
#    t = tfail
#    # getting desired debris catalog group
#    desiredCatalog = DR.desiredList(velcat,timecat,catalogList,Vmag,t)
#    # generating random debris pieces
#    massList,arefList,velList,ncdList,minfcdList,CDList,nclList,minfclList,CLList,arefMean,numberOfPiecesMeanList,blast=DR.generateDebris(desiredCatalog,numPiecesPerSample,debrisCatalogFilePATH)
#
##    minmass,maxmass,totalMass=debrisCatalog.getMassBounds(desiredCatalog)
#
#    #######################################################################################
#    ##########Debris Propagation ##########################################################
#
#    #        debrisNumTimeSteps = np.array([],dtype=np.int32)
#    debrisNumTimeSteps = []
#    debrisStorage = np.array([],dtype=np.float64)
#
#    debrisIndex = coeffIX
#
#    ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#    Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
#
#    # Putting these here in the hope that they don't fall out of scope after windIX == 0
#    stateVec = -1
#    thetag = -1
#    mass = -1
#    indexEvent = -1
#    Vmag = -1
#
##    # generating new atmos profile
##    densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
##    atmProfile1 = [altitudeList,densityList,uList,vList,wList]
##
##    # Generating new random state vector
##    stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile1,Toffset,ThrustOffsetAngDeg)
##    print 'stateVec = ' + str(stateVec)
#
#    # Generate all the atmosphere profiles and trajectories up front
#    atmStorage = []
#    stateVecStorage = []
#    thetagStorage = []    #This is a constant for a given time of launch and failure
#
#    if (numWindSamples == numTrajSamples):
#        for windIX in range(numWindSamples):
#            # generating new atmos profile
#            densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
#            atmStorage.append([altitudeList,densityList,uList,vList,wList])
#
#        for trajIX in range(numTrajSamples):
#            ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#            Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
##            stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmStorage[trajIX],Toffset,ThrustOffsetAngDeg)
#
#            stateVec,thetag = GenerateRandomStateVectorFullTJC(missionList, tLaunchDesired, dt, atmStorage[trajIX], Toffset, ThrustOffsetAngDeg)
#            singleStateVec = [stateVec[0][tfailIX], stateVec[1][tfailIX], stateVec[2][tfailIX], stateVec[3][tfailIX], stateVec[4][tfailIX], stateVec[5][tfailIX]]
#            stateVecStorage.append(singleStateVec)
#            thetagStorage.append(thetag[tfailIX])
##            print 'stateVec = ' + str(singleStateVec)
#
#    elif (numTrajSamples == 1):
#        for windIX in range(numWindSamples):
#            # generating new atmos profile
#            densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
#            atmStorage.append([altitudeList,densityList,uList,vList,wList])
#
#        ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#        Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
#                                                                                        # vv HERE vv #
#        stateVec,thetag = GenerateRandomStateVectorFullTJC(missionList, tLaunchDesired, dt, atmStorage[0], Toffset, ThrustOffsetAngDeg)
#        singleStateVec = [stateVec[0][tfailIX], stateVec[1][tfailIX], stateVec[2][tfailIX], stateVec[3][tfailIX], stateVec[4][tfailIX], stateVec[5][tfailIX]]
#        stateVecStorage.append(singleStateVec)
#        thetagStorage.append(thetag[tfailIX])
#
##        stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmStorage[0],Toffset,ThrustOffsetAngDeg)
##        stateVecStorage.append(stateVec)
##        thetagStorage.append(thetag)
#
#
#    elif (numWindSamples == 1):
#        # generating new atmos profile
#        densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
#        atmStorage.append([altitudeList,densityList,uList,vList,wList])
#
#        for trajIX in range(numTrajSamples):
#            ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#            Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
#
#            stateVec,thetag = GenerateRandomStateVectorFullTJC(missionList, tLaunchDesired, dt, atmStorage[0], Toffset, ThrustOffsetAngDeg)
#            singleStateVec = [stateVec[0][tfailIX], stateVec[1][tfailIX], stateVec[2][tfailIX], stateVec[3][tfailIX], stateVec[4][tfailIX], stateVec[5][tfailIX]]
#            stateVecStorage.append(singleStateVec)
#            thetagStorage.append(thetag[tfailIX])
#
##            stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmStorage[0],Toffset,ThrustOffsetAngDeg)
##            stateVecStorage.append(stateVec)
##            thetagStorage.append(thetag)
#
#    else:
#        print 'THIS OPTION IS NOT HANDLED YET a9ds8fy'
#        return
#
#
#
#
#
#
##    pdb.set_trace()
#
##    for windIX in range(numWindSamples):
#
#
##        # Don't generate a new wind profile on the first time through, use the one you've already got
##        if (windIX > 0):
##            densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
##            atmProfile1 = [altitudeList,densityList,uList,vList,wList]
##
##        for trajIX in range(numTrajSamples):
##
##            # Don't generate a new random state vector on the first time through
##            if (trajIX > 0) || (windIX > 0):
##                #random state vector... recalculated interpolation matrices every single time
##                stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile1,Toffset,ThrustOffsetAngDeg)
##                print 'stateVec = ' + str(stateVec)
#
#
##        # generating new atmos profile
##        densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
##        atmProfile1 = [altitudeList,densityList,uList,vList,wList]
##
##        if (windIX == 0):   #Only want to compute a single trajectory but have many wind profiles, so just do it for the first windIX
##            # random state vector... recalculated interpolation matrices every single
##            stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile1,Toffset,ThrustOffsetAngDeg)
##            print 'stateVec = ' + str(stateVec)
#
#    print 'Entering Debris Prop Area'
#
#    for windIX in range(numWindSamples):
#        altitudeList = atmStorage[windIX][0]
#        densityList = atmStorage[windIX][1]
#        uList = atmStorage[windIX][2]
#        vList = atmStorage[windIX][3]
#        wList = atmStorage[windIX][4]
#
#        for trajIX in range(numTrajSamples):
#            stateVec = stateVecStorage[trajIX]
#            thetag = thetagStorage[trajIX]
#
#
#        #    for debrisIndex in range(len(massList)):
#            ncdval = ncdList[debrisIndex]
#            minfcdval = minfcdList[debrisIndex]
#            nclval = nclList[debrisIndex]
#            minfclval = minfclList[debrisIndex]
#            CDval = CDList[debrisIndex]
#            CLval = CLList[debrisIndex]
#
#            for pieceIndex in range(numPiecesPerSample):
#
#                massval = massList[debrisIndex][pieceIndex]
#                arefval = arefList[debrisIndex][pieceIndex]
#                velImpval = velList[debrisIndex][pieceIndex]
#                #velImpval = 100.  # addition...just to check behavior
#                theta1 = np.random.uniform(0.0,2*np.pi)
#                theta2 = np.random.uniform(0.0,2*np.pi) #it is not required this to be 0 to 2pi, domain is covered with 0 to pi since prev angle helps cover the entire sphere
#                Vdeb = np.array([[velImpval,0,0]]).T # velocity impulse calculation
#                Rz = orbitTools.cRnMatrix_Zaxis(theta1)
#                Rx = orbitTools.cRnMatrix_Xaxis(theta2)
#                Vdeb = np.dot(Rz,Vdeb)
#                Vdeb = np.dot(Rx,Vdeb)
#                Vdeb = Vdeb[:,0] # making it a 1 by 3 vector...impulse velocity in cartesin coordinates
#                #print arefval
#                # debris propagation calculation
#
#    #            print 'pieceIX = ' + str(pieceIndex) + '   mass = ' + str(massval)
#
#                debrisResults, numFinalSteps = dp.debrispropagation(initialstate = stateVec,
#                                                                    debrisvel = Vdeb,
#                                                                    mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=1,
#                                                                    minfcl=minfclval,cl=CLval,loverd = 0,atmosoption=2,altitudelist=altitudeList,
#                                                                    densitylist=densityList,ulist=uList,vlist=vList,
#                                                                    wlist=wList,geoptions=0,filename='none',planetmodel=0,dtinterval = dt,ndtinterval = ndtinterval,thetag0=thetag,
#                                                                    ncd=ncdval,
#                                                                    ncl=nclval,
#                                                                    nlist=len(altitudeList))
#
##                pdb.set_trace()
#
#                altitudeFinal = debrisResults[6*numFinalSteps+2]
#                latitudeFinal = debrisResults[0]
#                longitudeFinal = debrisResults[1]
#                Vfinal = debrisResults[3] # final velocity magnitude relative to Earth
#                if altitudeFinal>0:
#                    print 'Warning...debris is not on the ground. Final altitude is ',altitudeFinal
#
#                #                pdb.set_trace()
#
#                # Saves the stuff I care about
#                #                debrisStorage.append(debrisResults[0:(numFinalSteps)*6].reshape(numFinalSteps,6)[:,:3])
#                # Chops off the tail of debrisResults, reshapes it into matrix form, then chops off the last three columns
#                # Then flattens out the result to concatentate it to the returning structure
#                # For some reason, those extra parentheses are necessary to avoid errors when concatenating the first (empty) array
#                debrisStorage = np.concatenate((debrisStorage,  (debrisResults[0:(numFinalSteps)*6].reshape(numFinalSteps,6)[:,:3]).flatten()) )
#                debrisNumTimeSteps.append(numFinalSteps)
#
#    # Convert the time steps array to be numpy int32 array
#    debrisNumTimeSteps = np.array(debrisNumTimeSteps,dtype=np.int32)
#    print 'DONE'
#
#    ###### Now to read these points into a PointCloud, for input into a footprint and export to GE
#    UTC = curMission['initialUTC'] + (tLaunchDesired + tfail)/(24*3600.)
#    zBinHeightKm = 5.0
#
#    launchLat = 28.445455
#    launchLon = -80.564865
#    launchAzimuth = 0.976192685683
#
#    all_points_delta_t = 1.0    #minutes
#
#    numRuns = len(debrisNumTimeSteps)
#    maxTime = debrisNumTimeSteps.max()*dt
#    myPointCloud = ceb.PyPointCloud(debrisStorage, numRuns, debrisNumTimeSteps, maxTime, curMission['deltaT'],
#                                    UTC, all_points_delta_t,
#                                    launchLat, launchLon, launchAzimuth)
#
#    # Pack things up in a dictionary (for pickling)
#    mpc = dict(flatPointArray = debrisStorage, numPieces = numRuns, numTimeSteps = debrisNumTimeSteps, maxTime = maxTime, deltaTsec = curMission['deltaT'],
#               UTC = UTC, all_points_delta_t = all_points_delta_t,
#               launchLat = launchLat, launchLon = launchLon, launchAzimuth = launchAzimuth)
#    return mpc
#    #return myPointCloud








# Just some general functions that I want hanging around in here
#(missionList,tLaunchDesired,tfail,dt,atmStorage[0],Toffset,ThrustOffsetAngDeg)
#def GenerateRandomStateVectorFullTJC(missionList, tLaunchDesired, dt, atmProfile, Toffset, ThrustOffsetAngDeg):
#
#    myTrajList = traj.generateRandomTrajFullFAST(missionList,tLaunchDesired,dt,atmProfile,Toffset,ThrustOffsetAngDeg)
#
#    myTrajStack = np.dstack(myTrajList)
#    myTrajStackFinal = []
#    for ix in range(9):
#        myTrajStackFinal.append(np.hstack( [myTrajStack[0][ix][0], myTrajStack[0][ix][1][1:], myTrajStack[0][ix][2][1:]]  ))
#
#    [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector] = myTrajStackFinal    #index 0=first stage, 1=coast period, 2=second stage
#
#    stateVector = [xc,yc,zc,Vxc,Vyc,Vzc]
#    return stateVector, thetagVector


#def GenerateRandomStateVectorFullTJC(missionList, tLaunchDesired, dt, atmProfile, Toffset, ThrustOffsetAngDeg):
#
#    myTrajList = traj.generateRandomTrajFullFAST(missionList,tLaunchDesired,dt,atmProfile,Toffset,ThrustOffsetAngDeg)
#
#    myTrajStack = np.dstack(myTrajList)
#    myTrajStackFinal = []
#    for ix in range(9):
#        myTrajStackFinal.append(np.hstack( [myTrajStack[0][ix][0], myTrajStack[0][ix][1][1:], myTrajStack[0][ix][2][1:]]  ))
#
#    [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector] = myTrajStackFinal    #index 0=first stage, 1=coast period, 2=second stage
#
#
#    stateVector = [xc,yc,zc,Vxc,Vyc,Vzc]
#    return stateVector, thetagVector





## Only does the first stage
#def GenerateRandomStateVectorFirstOnlyTJC(missionList, tLaunchDesired, dt, atmProfile, Toffset, ThrustOffsetAngDeg):
#
##    myTrajList = traj.generateRandomTrajFAST(missionList,tLaunchDesired,dt,atmProfile,Toffset,ThrustOffsetAngDeg)
#
#
#    myTrajList = traj.generateRandomTraj(missionList,tLaunchDesired,dt,atmProfile,Toffset,ThrustOffsetAngDeg)
#    [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent] = myTrajList[0]    #index 0=first stage, 1=coast period, 2=second stage
#
#    stateVector = [xc,yc,zc,Vxc,Vyc,Vzc]
#
#    return stateVector, thetagVector





##(atmoFile, trajectoryFiles, debrisCatalogFile, debrisCatalogFilePATH, dt, tLaunchDesired, tfail, nSamples, nRuns):
#def MonteCarloSingleBallisticCoeff__PROTOTYPE(curMission, coeffIX, tfail, numTrajSamples, numWindSamples, numPiecesPerSample):
#    print 'In MonteCarloSingleBallisticCoeff ======!!!!!PROTOTYPE!!!!=====, tLaunchDesired = ZERO'
#
#    # Unpack some of the mission variables
#    debrisCatalogFile = curMission['debrisCatPath'] + curMission['debrisCatFile']
#    debrisCatalogFilePATH = curMission['debrisCatPath']
#    dt = curMission['deltaT']
#
#    ndtinterval = 5*3600/dt     # is the upper bound for how long you think the debris propagation will run (5hrs)
#
#    # To read in
#    import pickle
#    [altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist] = pickle.load(open('mission1_atmo.pkl','rb'))
#
##    # This should eventually happen BEFORE this function or speed it up
##    atmoFile = curMission['atmosphere']
##    # getting atmospheric profile
##    altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints= AP.readGramAtmos(atmoFile)
##    # making sure atmospheric profile is decreasing in altitude
##    altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist= AP.maxAltitudeAtmosFix(altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints)
#
#
#    # These are some values that we're just going to hardcode
#    tLaunchDesired = 0.         # Delay in seconds of launch from nominal time
#    ndtinterval = 5*3600/dt     # is the upper bound for how long you think the debris propagation will run (5hrs)
#    Toffset_std = 0.01;         # Standard Dev of the thrust offset coeff (mean of 1)
#
#    # Get the mission list.  Eventually this will change format
#    missionList = loadmat(curMission['trajMatFile'])['missionList'] # loading trajectories  from SPOT
#
##    tfailIX = tfail/dt;
##    if (np.mod(tfail,dt) != 0.0):
##        print 'YOU PICKED A BAD FAILURE TIMESTEP'
##        return 0
##
#
#    #######debris catalog selection and generation#################################
#    # Need a trajectory first, which requires a wind profile, in order to get a Vmag.
#    # Make temporary ones of both and have them fall out of scope.  This is a hack.
#    # Hopefully future versions will have easily accessible Vmag's from the nominal solution
#    Vmag = -5.
#
#    if (1 == 1):
#        densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
#        atmProfile1 = [altitudeList,densityList,uList,vList,wList]
#        ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
#        Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
#        stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile1,Toffset,ThrustOffsetAngDeg)
#
#    # getting debris catalog information, getting desired debris catalog group
#    velcat,timecat,catalogList = DR.readDebris(debrisCatalogFile)
#    desiredCatalog = DR.desiredList(velcat,timecat,catalogList,Vmag,tfail)
#
#    # generating random debris pieces
#    massList,arefList,velList,ncdList,minfcdList,CDList,nclList,minfclList,CLList,arefMean,numberOfPiecesMeanList,blast=DR.generateDebris(desiredCatalog,numPiecesPerSample,debrisCatalogFilePATH)
#
#    #######################################################################################
#    ##########Debris Propagation ##########################################################
#
#    #        debrisNumTimeSteps = np.array([],dtype=np.int32)
#
#    # These two arrays will be appended to and returned as the solution of this function
#    debrisNumTimeSteps = []
#    debrisStorage = np.array([],dtype=np.float64)
#
#    # Only want to do a single debris coefficient, use the one passed in
#    if (coeffIX == []):
#        # If it's empty, then let's look at all the debris pieces
#        coeffIX = range(len(massList))
#    # If it's not empty, then just look at the pieces specified (probably -1)
##    debrisIndex = coeffIX
#
##    ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
##    Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
#
##    # Putting these here in the hope that they don't fall out of scope after windIX == 0
##    stateVec = -1
##    thetag = -1
##    mass = -1
##    indexEvent = -1
##    Vmag = -1
#
#    #    # generating new atmos profile
#    #    densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
#    #    atmProfile1 = [altitudeList,densityList,uList,vList,wList]
#    #
#    #    # Generating new random state vector
#    #    stateVec,thetag,mass,indexEvent,Vmag = traj.generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile1,Toffset,ThrustOffsetAngDeg)
#    #    print 'stateVec = ' + str(stateVec)
#
##    # Generate all the atmosphere profiles and trajectories up front
##    atmStorage = []
##    stateVecStorage = []
##    thetagStorage = []    #This is a constant for a given time of launch and failure
##    tfailStorage = []
##
##    # Generate the wind profiles first
##    for windIX in range(numWindSamples):
##        print 'windIX = ' + str(windIX)
##        # generating new atmos profile
##        densityList,uList,vList,wList= AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
##        atmStorage.append([altitudeList,densityList,uList,vList,wList])
##
##        stateVecStorage_TEMP = []
##        thetagStorage_TEMP = []
##        tfailStorage_TEMP = []
##        # For every wind profile, generate a trajectory
##        for trajIX in range(numTrajSamples):
##            ThrustOffsetAngDeg = [np.random.normal(0,.1),np.random.normal(0,.1)] # offset in thrust angles
##            Toffset = np.random.normal(1,Toffset_std) # value to multiply the thrust...
##
##            singleStateVec,thetag,tfail_real = \
##                GenerateRandomStateVectorNearestTimeFail(missionList, tLaunchDesired, dt, atmStorage[windIX], Toffset, ThrustOffsetAngDeg, tfail)
##
##            stateVecStorage_TEMP.append(singleStateVec)
##            thetagStorage_TEMP.append(thetag)
##            tfailStorage_TEMP.append(tfail_real)
##
##        stateVecStorage.append(stateVecStorage_TEMP)
##        thetagStorage.append(thetagStorage_TEMP)
##        tfailStorage.append(tfailStorage_TEMP)
#
#
##    return stateVecStorage
#
#
#    print 'Entering Debris Prop Area'
#
#    for windIX in range(numWindSamples):
#        altitudeList = atmStorage[windIX][0]
#        densityList = atmStorage[windIX][1]
#        uList = atmStorage[windIX][2]
#        vList = atmStorage[windIX][3]
#        wList = atmStorage[windIX][4]
#
#        for trajIX in range(numTrajSamples):
#            stateVec = stateVecStorage[windIX][trajIX]
#            thetag = thetagStorage[windIX][trajIX]
#
#
#            for debrisIndex in coeffIX:
#                ncdval = ncdList[debrisIndex]
#                minfcdval = minfcdList[debrisIndex]
#                nclval = nclList[debrisIndex]
#                minfclval = minfclList[debrisIndex]
#                CDval = CDList[debrisIndex]
#                CLval = CLList[debrisIndex]
#
#                for pieceIndex in range(numPiecesPerSample):
#
#                    massval = massList[debrisIndex][pieceIndex]
#                    arefval = arefList[debrisIndex][pieceIndex]
#                    velImpval = velList[debrisIndex][pieceIndex]
#                    #velImpval = 100.  # addition...just to check behavior
#                    theta1 = np.random.uniform(0.0,2*np.pi)
#                    theta2 = np.random.uniform(0.0,2*np.pi) #it is not required this to be 0 to 2pi, domain is covered with 0 to pi since prev angle helps cover the entire sphere
#                    Vdeb = np.array([[velImpval,0,0]]).T # velocity impulse calculation
#                    Rz = orbitTools.cRnMatrix_Zaxis(theta1)
#                    Rx = orbitTools.cRnMatrix_Xaxis(theta2)
#                    Vdeb = np.dot(Rz,Vdeb)
#                    Vdeb = np.dot(Rx,Vdeb)
#                    Vdeb = Vdeb[:,0] # making it a 1 by 3 vector...impulse velocity in cartesin coordinates
#                    #print arefval
#                    # debris propagation calculation
#
#                    #            print 'pieceIX = ' + str(pieceIndex) + '   mass = ' + str(massval)
#
#                    debrisResults, numFinalSteps = dp.debrispropagation(initialstate = stateVec,
#                                                                        debrisvel = Vdeb,
#                                                                        mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=1,
#                                                                        minfcl=minfclval,cl=CLval,loverd = 0,atmosoption=2,altitudelist=altitudeList,
#                                                                        densitylist=densityList,ulist=uList,vlist=vList,
#                                                                        wlist=wList,geoptions=0,filename='none',planetmodel=0,dtinterval = dt,ndtinterval = ndtinterval,thetag0=thetag,
#                                                                        ncd=ncdval,
#                                                                        ncl=nclval,
#                                                                        nlist=len(altitudeList))
#
#                    #                pdb.set_trace()
#
#                    altitudeFinal = debrisResults[6*numFinalSteps+2]
#    #                latitudeFinal = debrisResults[0]
#    #                longitudeFinal = debrisResults[1]
#                    Vfinal = debrisResults[3] # final velocity magnitude relative to Earth
#                    if altitudeFinal>0:
#                        print 'Warning...debris is not on the ground. Final altitude is ',altitudeFinal
#
#                    #                pdb.set_trace()
#
#                    # Saves the stuff I care about
#                    #                debrisStorage.append(debrisResults[0:(numFinalSteps)*6].reshape(numFinalSteps,6)[:,:3])
#                    # Chops off the tail of debrisResults, reshapes it into matrix form, then chops off the last three columns
#                    # Then flattens out the result to concatentate it to the returning structure
#                    # For some reason, those extra parentheses are necessary to avoid errors when concatenating the first (empty) array
#                    debrisStorage = np.concatenate((debrisStorage,  (debrisResults[0:(numFinalSteps)*6].reshape(numFinalSteps,6)[:,:3]).flatten()) )
#                    debrisNumTimeSteps.append(numFinalSteps)
#
#    # Convert the time steps array to be numpy int32 array
#    debrisNumTimeSteps = np.array(debrisNumTimeSteps,dtype=np.int32)
#    print 'DONE'
#
#    ###### Now to read these points into a PointCloud, for input into a footprint and export to GE
#    UTC = curMission['initialUTC'] + (tLaunchDesired + tfail)/(24*3600.)
#    zBinHeightKm = 5.0
#
#    numRuns = len(debrisNumTimeSteps)
#    maxTime = debrisNumTimeSteps.max()*dt
#    myPointCloud = ceb.PyPointCloud(debrisStorage, numRuns, debrisNumTimeSteps, maxTime, curMission['deltaT'],
#                                    UTC, curMission['all_points_delta_t'],
#                                    curMission['launchLat'], curMission['launchLon'], curMission['launchAzimuth'])
#
#    # Pack things up in a dictionary (for pickling)
#    mpc = dict(flatPointArray = debrisStorage, numPieces = numRuns, numTimeSteps = debrisNumTimeSteps, maxTime = maxTime, deltaTsec = curMission['deltaT'],
#               UTC = UTC, all_points_delta_t = curMission['all_points_delta_t'],
#               launchLat = curMission['launchLat'], launchLon = curMission['launchLon'], launchAzimuth = curMission['launchAzimuth'])
#    return mpc



#def GenerateRandomStateVectorNearestTimeFail(missionList, tLaunchDesired, dt, atmProfile, Toffset, ThrustOffsetAngDeg, tfail):
#    
#    myTrajList = traj.generateRandomTraj(missionList,tLaunchDesired,dt,atmProfile,Toffset,ThrustOffsetAngDeg)
#    
#    myTrajStack = np.dstack(myTrajList)
#    myTrajStackFinal = []
#    for ix in range(9):
#        myTrajStackFinal.append(np.hstack( [myTrajStack[0][ix][0], myTrajStack[0][ix][1][1:], myTrajStack[0][ix][2][1:]]  ))
#    
#    [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector] = myTrajStackFinal    #index 0=first stage, 1=coast period, 2=second stage
#    
#    nearestIX = (np.abs(tc-tfail)).argmin()
#    stateVector = [xc[nearestIX],yc[nearestIX],zc[nearestIX],Vxc[nearestIX],Vyc[nearestIX],Vzc[nearestIX]]
#    
#    return stateVector, thetagVector[nearestIX], tc[nearestIX]


## Just some general functions that I want hanging around in here
#def TestForSpeed(curMission, numRuns):
#    print 'In GenerateLRCH'
#
#    from scipy.io import loadmat
#    import orbitTools as ot
#
#    missionList = loadmat(curMission['trajMatFile'])['missionList'] # loading trajectories  from SPOT
#    trajStorage = np.array([])
#    trajNumTimeSteps = []
#    maxTime = 0;
#    tLaunchDesired = 0.0 # launch time in seconds from beginning of launch window
#
#    Toffset=1.0
#    ThrustOffsetAngDeg = [0.,0.]
#
#    for tx in range(numRuns):
#        print 'tx = ' + str(tx)
#        myTrajList = traj.generateRandomTrajFAST(missionList,tLaunchDesired,curMission['deltaT'],[],Toffset,ThrustOffsetAngDeg)
