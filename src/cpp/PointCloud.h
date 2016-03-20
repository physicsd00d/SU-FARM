//
//  PointCloud.h
//  Prop3Dof
//
//  Created by Thomas Colvin on 9/25/13.
//  Copyright (c) 2013 Thomas Colvin. All rights reserved.
//

#ifndef __Prop3Dof__PointCloud__
#define __Prop3Dof__PointCloud__

#define INTxx int

#include <cstdlib>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <map> // header file needed for to use MAP STL
using std::map;

#include <fstream>
using std::ofstream;
using std::ifstream;
using std::ios;

#include <string>
using std::string;

#include "Point.h"

//#define NASkm 18.289             // [km] This is the top of the NAS



class PointCloud{
private:

    
protected:
    void Destruct_all_points();

    // Timing info
	double all_points_delta_t;  //seconds
	int all_points_num_range;
    
    // Location info
    double all_points_launchLat, all_points_launchLon, all_points_launchAzimuth;
    
    // Discretization info
//    double xBinLength, yBinLength, zBinHeight;  //[km]

    ifstream allPointsFile;
    vector<int> num_all_points_at_tstep;	//This applies to the all_points vector, not to the footprint vector
    vector<Point> all_points_single; //[point index]

    // Use this to make sure loading at each timestep is used correctly
	int lastTimeStep;   //Check code to make sure this is still useful
    
    vector<vector<Point> > all_points_total;	//[timestep][point index]
    
    map<int, int> totalNumPointsPassedInPerID;  //totalNumPointsPassedInPerID[beta / ID] = numHere
    
    double NASkm;
    
//    void assemble_all_points_debris(void *flatPointArray, void *pointIdArray_in, void *massArray_in, void *areaArray_in,
//                                    int numPieces, void *numTimeSteps, int maxTime, double deltaT, double timeOffsetSec, double reactionTimeMinutes);
    void assemble_all_points_debris(vector<double> flatPointArray_in, vector<int> pointIdArray_in, vector<double> massArray_in,
                                    vector<double> areaArray_in, int numPieces, vector<int> numTimeSteps_in, int maxTime,
                                    double deltaTsec, double timeOffsetSec, double reactionTimeMinutes);


public:
    PointCloud();
    PointCloud(string pointsFile);
    PointCloud(PointCloud *myCloud);
//    PointCloud(void *flatPointArray, void *pointIdArray, int numPieces, void *numTimeSteps, int maxTime, double deltaT,
//               double all_points_UTC, double all_points_delta_t, double timeOffsetSec,
//               double all_points_launchLat_in, double all_points_launchLon_in, double all_points_launchAzimuth_in,
//               void *massArray, void *areaArray, double reactionTimeMinutes);
    PointCloud(vector<double> flatPointArray_in, vector<int> pointIdArray, int numPieces, vector<int> numTimeSteps_in, int maxTime, double deltaTsec,
               double all_points_UTC_in, double all_points_delta_t_in, double timeOffsetSec,
               double all_points_launchLat_in, double all_points_launchLon_in, double all_points_launchAzimuth_in,
               vector<double> massArray, vector<double> areaArray, double reactionTimeMinutes, double NASkm);
    ~PointCloud();
    
    void PrintAllPoints();

    
    double all_points_UTC;
    
    void setTimingInfo(double all_points_UTC, double all_points_delta_t, int all_points_num_range);
    void setLocationInfo(double all_points_launchLat, double all_points_launchLon, double all_points_launchAzimuth);

    int ChopAfterSeconds(int numSeconds);

    // Functions for loading points files
    int load_points_file_not_points(string allPointsFileName);
    void load_points_at_timestep(int tstep);
    int load_points_file_and_points(string allPointsFileName);

    // Return the points vectors
    int getNumPointsAtTstep(int tx);
    
    double getDeltaT();
    double getLaunchLat();
    double getLaunchLon();
    double getLaunchAzimuth();
//    int getNumRange();
    int getPointsRange();
    double getInitialUTC();
    
//    virtual int incorporateNewPoints(vector<vector<Point> > &total_points_at, double InitialUTC, double tstepMinutes, double launchLat, double launchLon, double launchAzimuth);
    
    virtual int identifyYourself();
    virtual double getZBinHeight();
    vector<vector<Point> > getAllPoints();
    map <int, int> getTotalNumPointsPassedInPerID();
    
//    void assemble_all_points_debris(void *thisthing);
//    void assemble_all_points_debris_double(void **thisthing);

    void ExportPointsToMatlab(string fileName, double deltaZkm);




};



#endif /* defined(__Prop3Dof__PointCloud__) */
