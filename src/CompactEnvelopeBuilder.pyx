from libcpp.string cimport string
from libcpp cimport bool
import numpy as np
cimport numpy as np
#from ctypes import *
import ctypes
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.utility cimport pair

from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

#  http://stackoverflow.com/questions/19057439/missing-numpy-arrayobject-h-while-compiling-pyx-file
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC

from cpython.buffer cimport Py_buffer
from libc.string cimport memcpy

import datetime as dt


#cdef extern from "classes.h":
#    cdef cppclass BaseClass:
#        BaseClass(BaseClass *parent)
#    
#    cdef cppclass MainClass:
#        MainClass(BaseClass *parent)
#        void accumulate(int new_number)
#        int get_total()
#
#
#
### ~~~~~~~~~~~~~~~~~~~~~ SKYGRID CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#cdef extern from "PointCloud.h" namespace "PointCloud":
#    




## ~~~~~~~~~~~~~~~~~~~~~ SKYGRID CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cdef extern from "SkyGrid.h" namespace "SkyGrid":
cdef extern from "SkyGrid.h":
    cdef cppclass SkyGrid:
#        SkyGrid(char *CapeLrhcFile, double xBinLength_in, double yBinLength_in, double zBinHeight_in)
        SkyGrid(PointCloud *newCloud, double xBinLength_in, double yBinLength_in, double zBinHeight_in)

        double getLaunchLat()
        double getLaunchLon()
        double getLaunchAzimuth()
        double getInitialUTC()

#        void PythonDebrisIntoGrid(void *flatPointArray, int numPieces, void *numTimeSteps, int maxTime, double deltaT,  # arguments for assembling the points
#                                  double NewInitialUTC, double timeOffsetSec, double launchLat, double launchLon, double launchAzimuth)  #encorporation
        
        void PythonDebrisIntoGrid(PointCloud *incomingCloud)


#        void ExportBinnedDebrisGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min)

#        void generateAllPointsFromKDE(double deltaXY)
        
        void generateAllPointsFromGrid()
        
        int identifyYourself()
        int getNumRange()

        
        # Probability stuff
#        void ConvertToEmptyProbability()
        void ConvertToProbability(double weight)
#        void weightedCombine(SkyGrid *newSkyGrid, double weight)
#        bool isTotalProbabilityGood()

#        double generateAllPointsFromSimpleHistogram(double thresh, int Ntotal, int numEventsSimulated,  double pFail)

#        double generateAllPointsFromASH(vector[int] numberOfPiecesMeanArray_in, vector[double] arefMeanList_in, int numDebIX, double thresh, double pFail)
        
        
        void ASH2(double h1, double h2)
        double WTF()
        void DumpGridToMatlab(char *fileName)

        vector[double] createEmptyAircraftDensityMap()
        void populateAircraftDensityMap(void *densityMapArray, int numElements)

        void UploadAircraftTrackMap(map[int,pair[vector[vector[double]], string]] AircraftTrackMap)
        void UploadAircraftPropertiesMap(map[string,map[string,double]] AircraftPropertiesMap)

        map[int, double] CalculateRiskToIndividualAircraft(vector[int] numberOfPiecesMeanList, vector[double] arefMeanList, int secondsFromMidnightUTC)
        map[int, double] CalculateRiskToIndividualAircraft_OnTheFly(vector[int] numberOfPiecesMean, vector[double] arefMean, int secondsFromMidnightUTC,
                                                                     double h1_in, double h2_in)
        
#        vector[vector[double]] SendGridToPython(int tx_desired)
        map[double, map[double, map[double,double]]] SendGridToPython(int tx_desired)

        # For debugging
        map[int, map[int, map[int,double]]] getSpatialProbabilty()
        void generateSpatialProbability(int whichProb)
        map[int, map[int, map[int,double]]] projectSpatialProbabilityFAA(double newDeltaXY, double newDeltaZ)

        void generateHazardProbabilities(vector[int] numberOfPiecesMean)
        # double generateAllPoints_CumulativeTJC(double thresh, int whichProb)
        double generateAllPoints_CumulativeFAA(double thresh, int whichProb, double pFail)




#        np.ndarray[np.double_t, ndim=2] createEmptyAircraftDensityMap()





cdef class PySkyGrid:
    cdef SkyGrid *thisptr                    # hold a C++ instance which we're wrapping
    
    def __cinit__(self, PyPointCloud incoming, double xBinLength_in, double yBinLength_in, double zBinHeight_in):
        self.thisptr = new SkyGrid( incoming.thisptr,  xBinLength_in,  yBinLength_in,  zBinHeight_in)
    
    def __dealloc__(self):
        del self.thisptr

    def getLaunchLat(self):
        return self.thisptr.getLaunchLat()

    def getLaunchLon(self):
        return self.thisptr.getLaunchLon()
    
    def getLaunchAzimuth(self):
        return self.thisptr.getLaunchAzimuth()

    def getInitialUTC(self):
        return self.thisptr.getInitialUTC()
    
    def getNumRange(self):
        return self.thisptr.getNumRange()
    
    def IncorporatePointCloudIntoGrid(self, PyPointCloud incomingCloud):
        self.thisptr.PythonDebrisIntoGrid(incomingCloud.thisptr)

#    def ExportToGoogleEarth(self, char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min):
#        self.thisptr.ExportBinnedDebrisGoogleEarth(googleEarthFile, yyyy, mm, dd, hour, min)

#    def generateAllPointsFromKDE(self, double deltaXY):
#        self.thisptr.generateAllPointsFromKDE(deltaXY)

#    def generateAllPointsFromGrid(self):
#        self.thisptr.generateAllPointsFromGrid()
    
    def IdentifyYourself(self):
        self.thisptr.identifyYourself()

#    def ConvertToEmptyProbability(self):
#        self.thisptr.ConvertToEmptyProbability()

#    def ConvertToProbability(self,  double weight):
#        self.thisptr.ConvertToProbability(weight)

#    def weightedCombine(self, PySkyGrid newSkyGrid, double weight):
#        self.thisptr.weightedCombine(newSkyGrid.thisptr, weight)

#    def isTotalProbabilityGood(self):
#        return self.thisptr.isTotalProbabilityGood()

#    def generateAllPointsFromSimpleHistogram(self, double thresh, int Ntotal, int numEventsSimulated,  double pFail):
#        return self.thisptr.generateAllPointsFromSimpleHistogram(thresh, Ntotal, numEventsSimulated, pFail)
    
    def generateASH(self, h1, h2):
        self.thisptr.ASH2(h1, h2)

    
#    def generateAllPoints(self, numberOfPiecesMeanList, arefMeanList, double thresh, double pFail):
#        return self.thisptr.generateAllPointsFromASH(numberOfPiecesMeanList, arefMeanList, len(numberOfPiecesMeanList), thresh, pFail)
    
#    # Doesn't look like this gets used anymore (we call ASH separately at the moment.  Comment out.
#    def generateAllPointsFromASH(self, numberOfPiecesMeanList, arefMeanList, double thresh, double pFail,
#                                 double h1, double h2, int numDebrisPerIXSimulated):
#
#        self.thisptr.ASH2(h1, h2, numDebrisPerIXSimulated)
#        return self.thisptr.generateAllPointsFromASH(numberOfPiecesMeanList, arefMeanList, len(numberOfPiecesMeanList), thresh, pFail)

    def DumpGridToMatlab(self, char *filename):
        self.thisptr.DumpGridToMatlab(filename)

    # TODO: If the code is working without this stuff, then send it to a graveyard.
#     def createEmptyAircraftDensityMap(self):
# #        cdef vector<double> tempAns = self.thisptr.createEmptyAircraftDensityMap()

#         # i cannot figure out how to do this nicely, so it's going to become a slow-ass hack
#         tempAns = self.thisptr.createEmptyAircraftDensityMap()
#         numEntries = tempAns.size()
        
#         # Allocate the memory
#         latlonArray = np.zeros( (numEntries,))
#         for ix in range(numEntries):
#             latlonArray[ix] = tempAns[ix]

#         latlonArray = latlonArray.reshape(numEntries/4, 4)
#         return latlonArray

#     def populateAircraftDensityMap(self, np.ndarray[double, ndim = 1, mode="c"] densityMapArray, int numElements):
#         if (numElements > 0) or (numElements == -1):
#             self.thisptr.populateAircraftDensityMap(&densityMapArray[0], numElements)

    def SendGridToPython(self, int tx_desired):
        # Get the grid
        # Explicitly let cython know that the outer layer is going to be a dict/map and not a pair
        probGrid = <dict> self.thisptr.SendGridToPython(tx_desired)
                
        tempGrid = []
        for zval in probGrid:
            for lonVal in probGrid[zval]:
                for latVal in probGrid[zval][lonVal]:
                    probVal = probGrid[zval][lonVal][latVal]
                    tempGrid.append([zval, lonVal, latVal, probVal])
        probGrid = np.array(tempGrid)
        return probGrid

        
#        # You have to tell Cython what to expect for the sub-iterators (so also defining curz for consistency)
#        cdef map[double, map[double, map[double,double]]].iterator curz
#        cdef map[double, map[double,double]].iterator curx
#        cdef map[double,double].iterator cury
#
#        probGrid = []                                           # This is the storage list which we will return
#        
#        curz = curGrid.begin()                                  # iterator over z-values
#        while curz != curGrid.end():                            # loop through iterators until you hit the end
#             
#            z_alt = float(deref(curz).first)                    # first is the altitude in km
#            curx = (deref(curz).second).begin()                 # iterator over x-values
#            while curx != (deref(curz).second).end():
#                
#                longitude = float(deref(curx).first)            # first is longitude in degrees
#                cury = deref(curx).second.begin()               # iterator over y-values
#                while cury != deref(curx).second.end():
#                    
#                    latitude = float(deref(cury).first)         # first is latitude in degrees
#                    probHere = float(deref(cury).second)        # second is probability of strike
#                    
#                    probGrid.append([z_alt, longitude, latitude, probHere])
#                
#                    inc(cury)                                   # all of these increments are PRE increments (see imports at top)
#                inc(curx)
#            inc(curz)
#        return probGrid

    def UploadAircraft(self, ACdict):
#        # Just take one aircraft for now
#        keys = ACdict.keys()
#        singleAC = ACdict[keys[0]]
#        
#        testList = range(9)
#        cdef vector[double] ACvec = testList
#
#        cdef map[int,vector[vector[double]]] testMap = ACdict
#        self.thisptr.UploadSingleAircraft(testMap)
        self.thisptr.UploadAircraftTrackMap(ACdict)

    def UploadAircraftPropertiesMap(self, incomingMap):
        self.thisptr.UploadAircraftPropertiesMap(incomingMap)
        
        
    def CalculateRiskToIndividualAircraft(self, numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC):
        return self.thisptr.CalculateRiskToIndividualAircraft(numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC)

    def CalculateRiskToIndividualAircraft_OnTheFly(self, numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC,
                                                   h1_in, h2_in):
        return self.thisptr.CalculateRiskToIndividualAircraft_OnTheFly(numberOfPiecesMeanList, arefMeanList, secondsFromMidnightUTC,
                                                                       h1_in, h2_in)


    def GenerateSpatialProbability(self, whichProb):
        #define PROB_IMPACT      1001
        #define PROB_CASUALTY    1002
        #define PROB_CATASTROPHE 1003
        self.thisptr.generateSpatialProbability(whichProb)
            
    def GetSpatialProbabilty(self):
        return self.thisptr.getSpatialProbabilty()

    def GetSpatialProbabilty_Coarse(self, newDeltaXY, newDeltaZ):
        return self.thisptr.projectSpatialProbabilityFAA(newDeltaXY, newDeltaZ)

    def generateHazardProbabilities(self, numberOfPiecesMean):
        self.thisptr.generateHazardProbabilities(numberOfPiecesMean)

    # def generateAllPoints_CumulativeTJC(self, double thresh, int whichProb):
    #     return self.thisptr.generateAllPoints_CumulativeTJC(thresh, whichProb)
            
    def generateAllPoints_CumulativeFAA(self, double thresh, int whichProb, double pFail):
        return self.thisptr.generateAllPoints_CumulativeFAA(thresh, whichProb, pFail)



# This is here IN THE HOPES that it will allow me to pass PointCloud / SkyGrid objects into the footprint class
# http://docs.cython.org/src/userguide/external_C_code.html

cdef extern from "PointCloud.h":
# cdef extern from "PointCloud.h" namespace "PointCloud":
    cdef cppclass PointCloud:
        PointCloud()
#        PointCloud(void *flatPointArray_in, void *pointIdArray, int numPieces, void *numTimeSteps_in, int maxTime, double deltaTsec,
#                   double all_points_UTC, double all_points_delta_t, double timeOffsetSec,
#                   double all_points_launchLat_in, double all_points_launchLon_in, double all_points_launchAzimuth_in,
#                   void *massArray, void *areaArray, double reactionTimeMinutes)

        PointCloud(vector[double] flatPointArray_in, vector[int] pointIdArray, int numPieces, vector[int] numTimeSteps_in, int maxTime, double deltaTsec,
                   double all_points_UTC_in, double all_points_delta_t_in, double timeOffsetSec,
                   double all_points_launchLat_in, double all_points_launchLon_in, double all_points_launchAzimuth_in,
                   vector[double] massArray, vector[double] areaArray, double reactionTimeMinutes, double NASkm)
    
    
        void ExportPointsToMatlab(char *fileName, double deltaZkm)
        int ChopAfterSeconds(int numSeconds)


cdef class PyPointCloud:
    cdef PointCloud *thisptr                    # hold a C++ instance which we're wrapping
    
    def __cinit__(self):
        self.thisptr = new PointCloud()
    
    def __cinit__(self, dict pcd, double secondsFromLaunch, dict curMission):
        
        reactionTimeMinutes     = curMission['reactionTimeMinutes']
        all_points_delta_t      = curMission['all_points_delta_t']
        NASkm                   = curMission['NASkm']
        
        self.thisptr = new PointCloud(pcd['flatPointArray'],    pcd['debrisID'],        pcd['numPieces'],       pcd['numTimeSteps'],
                                      pcd['maxTime'],           pcd['deltaTsec'],
                                      pcd['UTC'],               all_points_delta_t,     secondsFromLaunch,
                                      pcd['launchLat'],         pcd['launchLon'],       pcd['launchAzimuth'],
                                      pcd['debrisMass'],        pcd['debrisArea'],      reactionTimeMinutes,    NASkm)
    
    def __dealloc__(self):
        del self.thisptr

    def ExportPointsToMatlab(self, char *fileName, double deltaZkm):
        self.thisptr.ExportPointsToMatlab(fileName, deltaZkm)

    def ChopAfterSeconds(self, int numSeconds):
        return self.thisptr.ChopAfterSeconds(numSeconds)



## ~~~~~~~~~~~~~~~~~~~~~ FOOTPRINT CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cdef extern from "Footprint3D.h":
    cdef cppclass Footprint3D:
        #        Footprint3D()
        Footprint3D(char *FootprintVectorFile)
        #Footprint3D(char *pointsFile, double bin_size_in)
        #void generate_footprint_at_timesteps()
#        Footprint3D(PointCloud &incoming)
        Footprint3D(PointCloud *incoming, double bin_size_in)
        
        void store_footprint_as_vector(char *footprintFileName)
        int load_footprint_as_vector(char *footprintFileName)
        
        void setAzimuthDeg(double azimuth_in)
        void setLaunchLatDeg(double lat)
        void setLaunchLonDeg(double lon)
        
        double getLaunchLat()
        double getLaunchLon()
        double getLaunchAzimuth()
        
        double getNumRange()
        double getDeltaT()
        
        void ChangeAzimuthToDeg(double newAzimuth)
        void ChangeLaunchSiteToDeg(double gdlatIN, double lonIN)
        void ShiftFootprintByMinutes(int HowManyMinutes)

        void exportGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min)
        void make_facet_files(char *folderName, int startTimeMinutes, int offsetTimeMinutes, int tstepMinutes)

        void MergeFootprintVectors(Footprint3D *incomingFP)
        void SmoothedOut(double newDeltaT)
        int ProjectAllPointsDown()
        double ChopTimeAt(double seconds)

        void SlideFootprintBySeconds(int howManySeconds)


#        void exportGoogleEarth(char *googleEarthFile)




cdef class PyFootprint:
    cdef Footprint3D *thisptr                    # hold a C++ instance which we're wrapping
    
    # Typeless overloaded constructor
    def __cinit__(self, incoming, asVector=False):
        if (asVector == False):
            bin_size_in = -5    # This is for a pure PointCloud, which is kind of no longer allowed.  Must be SkyGrid.
            # Must first cast the incoming python object to be a known python object PySkyGrid, then can cast its thisptr to PointCloud
            self.thisptr = new Footprint3D(<PointCloud*> (<PySkyGrid>incoming).thisptr, bin_size_in)
        else:
            self.thisptr = new Footprint3D(incoming)    #incoming is a string here
    
    def __dealloc__(self):
        del self.thisptr
    
    def StoreFootprintAsVector(self, char *footprintFileName):
        self.thisptr.store_footprint_as_vector(footprintFileName)
        
    def LoadFootprintAsVector(self, char *footprintFileName):
        return self.thisptr.load_footprint_as_vector(footprintFileName)
    
    def ExportGoogleEarth(self, char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min):
        self.thisptr.exportGoogleEarth(googleEarthFile, yyyy, mm, dd, hour, min)

    def ChangeAzimuthToDeg(self, double newAzimuth):
        self.thisptr.ChangeAzimuthToDeg(newAzimuth)

    def ChangeLaunchSiteToDeg(self, double newLat, double newLon):
        self.thisptr.ChangeLaunchSiteToDeg(newLat,newLon)
    
    def SetAzimuthDeg(self, double newAzimuth):
        self.thisptr.setAzimuthDeg(newAzimuth)
    
    def SetLaunchLatDeg(self, double newLat):
        self.thisptr.setLaunchLatDeg(newLat)
    
    def SetLaunchLonDeg(self, double newLon):
        self.thisptr.setLaunchLonDeg(newLon)
                
    def GetAzimuthDeg(self):
        return self.thisptr.getLaunchAzimuth()
    
    def GetLaunchLatDeg(self):
        return self.thisptr.getLaunchLat()
    
    def GetLaunchLonDeg(self):
        return self.thisptr.getLaunchLon()
    
    def getNumRange(self):
        return self.thisptr.getNumRange()
    
    def getDeltaT(self):
        return self.thisptr.getDeltaT()

    def ShiftFootprintByMinutes(self, int HowManyMinutes):
        self.thisptr.ShiftFootprintByMinutes(HowManyMinutes)
    
    def MergeFootprintVectors(self, PyFootprint incomingFP):
        self.thisptr.MergeFootprintVectors(incomingFP.thisptr)
    
    def SmoothedOut(self, newDeltaT = -1.):
        self.thisptr.SmoothedOut(newDeltaT)
    
    def ProjectAllPointsDown(self):
        return self.thisptr.ProjectAllPointsDown()
    
    def ChopTimeAt(self, double seconds):
        return self.thisptr.ChopTimeAt(seconds)
        
    def MakeFacetFiles(self, folderName, int startTimeMinutes, int offsetTimeMinutes, int tstepMinutes):
        print folderName
        self.thisptr.make_facet_files(folderName, startTimeMinutes, offsetTimeMinutes, tstepMinutes)

    def SlideFootprintBySeconds(self, howManySeconds):
        self.thisptr.SlideFootprintBySeconds(howManySeconds)




## ~~~~~~~~~~~~~~~~~~~~~ TRAJECTORY CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This class exists here solely to turn a time/lon/lat/alt file into a google earth animation
#
# This class used to get used for all of the trajectory propagations, but now I just use
#   Francisco's fortran file and do all of the controlling from python/cython.
#
#
cdef extern from "Trajectory.h":
    cdef cppclass Trajectory:
        Trajectory()
        
        # These files read in the trajectories to plot
        void loadPrecomputedFile(string inFileName, bool isDegrees, bool isAltMeters)
        void loadDebris(vector[double] flatPointArray, vector[int] pointIdArray, vector[int] numTimeSteps, int maxTime, double deltaTsec, double timeOffsetSec, double reactionTimeMinutes)

        void setLaunchTime(int launch_year_in, int launch_month_in, int launch_day_in, int launch_hours_in, int launch_minutes_in, int launch_seconds_in)
        int write_to_google_earth_native(string basename, int printThisMany)



cdef class PyTrajectory:
    cdef Trajectory *thisptr                    # hold a C++ instance which we're wrapping

    # This takes a file of the format time/lon/lat/alt
    def __cinit__(self):
        self.thisptr = new Trajectory()
    
    def loadPrecomputedTrajectory(self, string inFileName, bool isDegrees, bool isAltMeters):
        self.thisptr.loadPrecomputedFile(inFileName, isDegrees, isAltMeters)
    
    def loadDebrisTrajectory(self, dict pcd, double secondsFromLaunch, dict curMission):
        reactionTimeMinutes     = curMission['reactionTimeMinutes']
        self.thisptr.loadDebris(pcd['flatPointArray'],    pcd['debrisID'],      pcd['numTimeSteps'],
                                pcd['maxTime'],           pcd['deltaTsec'],
                                secondsFromLaunch,        reactionTimeMinutes)
    

    def ExportGoogleEarth(self, outFileName, curDateTime):
        if curDateTime == []:
            # Use some default values that won't be confused for real dates
            yyyy = 1984
            mon  = 5
            day  = 27
            hr   = 12
            min  = 30
            sec  = 00
        else:
            yyyy = curDateTime.year
            mon  = curDateTime.month
            day  = curDateTime.day
            hr   = curDateTime.hour
            min  = curDateTime.minute
            sec  = curDateTime.second

        self.thisptr.setLaunchTime(yyyy, mon, day, hr, min, sec)
        self.thisptr.write_to_google_earth_native(outFileName, 1)






## ~~~~~~~~~~~~~~~~~~~~~ POINT CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This class exists solely for debugging purposes.
# You should never actually need to use this directly from Python.
#
#
cdef extern from "Point.h":
    cdef cppclass Point:
        Point(double gdLatIN, double lonIN, double z_in, double R_local_in)    #(rad, rad, km, km)

cdef class PyPoint:
    cdef Point *thisptr                    # hold a C++ instance which we're wrapping
    
    # Constructor.  R_local is no longer used; use anything.
    def __cinit__(self, gdLatIN, lonIN, z_in):
        R_local = -1
        self.thisptr = new Point(gdLatIN, lonIN, z_in, R_local)

















    