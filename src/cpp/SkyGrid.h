//
//  SkyGrid.h
//  Prop3Dof
//
//  Created by Thomas Colvin on 9/25/13.
//  Copyright (c) 2013 Thomas Colvin. All rights reserved.
//

#ifndef __Prop3Dof__SkyGrid__
#define __Prop3Dof__SkyGrid__

#include <iostream>

#include <stdio.h>


//#include <vector>
//using std::vector;

//#include "Point.h"
#include "Footprint3D.h"    //Contains a lot of the includes that this needs
#include "PointCloud.h"

//#include <gsl/gsl_blas.h>   //WHY???

#include <map> // header file needed for to use MAP STL
using std::map;

#include <utility>
using std::pair;

#include <algorithm>    // std::sort
using std::max;

// Options are PROB_IMPACT, PROB_CASUALTY, PROB_CATASTROPHE
#define PROB_IMPACT      1001
#define PROB_CASUALTY    1002
#define PROB_CATASTROPHE 1003

#define STORE_IX -666


#define XREF 0.
#define YREF 0.
#define ZREF 0.

class SkyGrid: public PointCloud{
private:
    double xBinLength, yBinLength, zBinHeight;  //[km]
    
    
    // The structure that keeps track of the quantities of interest in every bin cell
    struct binData{
        double probDebris;   // If used in GridMapDebIX this is effectively an int, counts pieces
        // If used in ProbabilityMapDebIX then it is the probability
        double probImpact;
        double probCasualty;
        double probCatastrophe;
        
        // Transitioning to these.  Delete the above eventually when they're fully phased out.
        double probNoImpact;
        double probNoCasualty;
        double probNoCatastrophe;
        
        double minVel;  // Velocity stats are measured in km/s
        double maxVel;
        double avgVel;
        
        double minMass;  // Mass stats are measured in kg
        double maxMass;
        double avgMass;
        
        double minArea;  // Area stats are measured in km^2
        double maxArea;
        double avgArea;
        
        int numTouched; // Used for spreading averages around.  Keeps track of how many times a cell has been added to during ASH.
        
        binData() {
            probDebris = 0;
            probImpact = 0;
            probCasualty = 0.;
            probCatastrophe = 0.;

            probNoImpact = 0;
            probNoCasualty = 0.;
            probNoCatastrophe = 0.;
            
            minVel = 1e10;
            maxVel = 0;
            avgVel = 0;
            
            minMass = 1e10;  // Mass stats are measured in kg
            maxMass = 0;
            avgMass = 0;
            
            minArea = 1e10;  // Area stats are measured in km^2
            maxArea = 0;
            avgArea = 0;
            
            numTouched = 0;
        }
        
    };
    
    map<int, int> numPtsLessThanReactionInNAS;  // Keep track of this because it will be needed for checking that the probabilities are consistent

    map<int, map<int, map<int, double> > > ACDensityMap;

    map<int, map<int, map<int, map<int, map<int,binData> > > > >GridMapDebIX;            //[t][debIX][z][x][y]   NOT THE SAME
    map<int, map<int, map<int, map<int, map<int,binData> > > > >ProbabilityMapDebIX;     //[t][z][x][y][debIX]   NOT THE SAME
    map<int, map<int, int> >totalNumPtsAtStepMapDebIX;   //[t][debIX]
    
    // This is a map that contains the probabilities of interest (as specified by whichProb) that have been
    //  summed over all time and over all debIX.  Thus, it represents that total probability of, for instance, impact
    //  at a given point in xyz space.
    map<int, map<int, map<int,double> > > SpatialProbabilty;

    vector<double> xpts_KDE;
    vector<double> ypts_KDE;
    
//    vector<vector<vector<vector<double> > > > Probability;
//    map<int, map<int, map<int, map<int,double> > > > ProbabilityMap;


    double hist_coeff;                      //multiply points in bin by this coeff to get probability
    
    // Some constants of the architecture
//    const static double NASkm = 18.289;             // [km] This is the top of the NAS
    const static double bufferDist = 0.;          // Want at least 10km on any edge for a histogram buffer.  TURNING OFF BUFFER!!!!
    int halfBufferCells;                                    // The number of extra cells you need on a side to have the bufferDist
    
    bool isProbability;     // Once you convert the object to a probability, you can never go back and add more grid points.
        
//    void GridTheSky(vector<vector<Point> > &total_points_at);   //[tx][point]   Only gets called in the constructorss
    void GridTheSky();   //[tx][point]   Only gets called in the constructorss

//    vector<vector<double> > generateProbabilityOfImpact(int tx, vector<int> numberOfPiecesMean, double thresh, double pFail);
    
    void loadRemainingPointsIntoAllPoints(int tx, vector<vector<double> > ProbabilityHere);


    
    map<int, pair<vector<vector<double> >, string> > AircraftTrackMap;
    map<string,map<string,double> > AircraftPropertiesMap;

    
    double generateAllPointsFromProbability(vector<int> numberOfPiecesMean, vector<double> arefMeanList, double thresh, double pFail);
//    double generateAllPointsFromProbability(double thresh, int Ntotal, int numEventsSimulated,  double pFail);
//    double generateAllPointsFromProbability002(double thresh, int Ntotal, int numEventsSimulated, double pFail);

    bool doneASH;
    bool hazardProbabilitiesGenerated;

    
public:
    // ==== Constructors (in order of relevance, only first one ever gets used anymore) ====
    SkyGrid(PointCloud *newCloud, double xBinLength_in, double yBinLength_in, double zBinHeight_in);
    SkyGrid(string CapeLrhcFile, double xBinLength_in, double yBinLength_in, double zBinHeight_in);
    SkyGrid(string SkyGridFileName);
    

    
    void RemoveHistogramKeepTiming();

        
    void ExportBinnedDebrisGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min);
    void StoreGrid(string outFileName);
    
    void ASH2(double h1, double h2);
    
    void UploadAircraftTrackMap(map<int, pair<vector<vector<double> >, string> > AircraftTrackMap_in, int aircraftTrackDeltaTSec);
//    void UploadAircraftTrackMap(map<int, vector<vector<double> > > AircraftTrackMap);
    void UploadAircraftPropertiesMap(map<string,map<string,double> > AircraftPropertiesMap_in);
    // map<int, double> CalculateRiskToIndividualAircraft(vector<int> numberOfPiecesMeanList, vector<double> arefMeanList, int secondsFromMidnightUTC);

    map<int, double> CalculateRiskToIndividualAircraft_OnTheFly(vector<int> numberOfPiecesMeanList, vector<double> arefMeanList, int secondsFromMidnightUTC,
                                                                double h1_in, double h2_in);
    vector<map<int,binData> > ASHDesiredPoint(double h1_in, double h2_in, vector<vector<double> > desiredPt);

    vector<double> ProbNoConsequence(map<int, binData> &probBeta, vector<int> numberOfPiecesMean,
                                     double cellVolume, double d_Airplane_top,
                                     double d_Airplane_front, double aircraftSpeed, double delta_t);
    
    // Debugging Functions (Matlab)
    void GoMatlab(string fileName, vector<Point> tempVec);
    void DumpGridToMatlab(char *fileName);
    
    // Debugging Functions (Used for being interactive with Python)
    map<int, map<int, map<int,double> > > getSpatialProbabilty();
    void generateSpatialProbability(int whichProb, int J_maxTimeStep, int f_startTimeStep);
    map<int, map<int, map<int,double> > > projectSpatialProbabilityFAA(double newDeltaXY, double newDeltaZ);

    
    
    void ConvertToEmptyProbability();
    void ConvertToProbability(double weight);
    void weightedCombine(SkyGrid *newSkyGrid, double weight);
    bool isTotalProbabilityGood();

    
    vector<Point> getGridValues(int tstep, int zstep);
    void generateAllPointsFromGrid();
    
    double generateAllPointsFromSimpleHistogram(double thresh, int Ntotal, int numEventsSimulated,  double pFail);
//    double generateAllPointsFromASH(vector<int> numberOfPiecesMeanArray_in, vector<double> arefMeanList_in, int numDebIX, double thresh, double pFail);
    

    void generateHazardProbabilities(vector<int> numberOfPiecesMean);
    double generateAllPoints_CumulativeFAA(double thresh, int whichProb, double pFail);
    
    


    map<int, vector<vector<double> > > ProbabilityTotalStorage;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
    map<int, int> stopIXStorage;
    
//    vector<vector<double> >SendGridToPython(int tx_desired);
    map<double, map<double, map<double,double> > >SendGridToPython(int tx_desired);
    
    int identifyYourself();
    double getZBinHeight();
    int getNumRange();
//    bool get_isProbability();

    // This function gets called from Python, assumes an existing grid is in place, and handles everything from there
    void PythonDebrisIntoGrid(void *flatPointArray, int numPieces, void *numTimeSteps, int maxTime, double deltaT,    // arguments for assembling the points
                              double NewInitialUTC, double timeOffsetSec, double launchLat, double launchLon, double launchAzimuth);    // and encorporating into grid
    void PythonDebrisIntoGrid(PointCloud *incomingCloud);


};


#endif /* defined(__Prop3Dof__SkyGrid__) */
