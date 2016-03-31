/*
 *  Footprint3D.h
 *  CompactEnvelope3D
 *
 *  Created by Thomas Colvin on 3/20/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FOOTPRINT_CLASS
#define FOOTPRINT_CLASS

#define BIN2D		1
#define PLAIN	2

#define DU 6378.145		// km  
#define R_equator 6378.1370   //[km]
#define R_polar 6356.7523	  //[km]

#include <cstring>
//#include <string>
//using std::string;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::ofstream;
using std::ifstream;
using std::ios;

#include <iomanip>

#include <list>
using std::list;

#include <vector>
using std::vector;

#include "Point.h"
#include "PointCloud.h"
#include "SkyGrid.h"


#include <cmath>
using std::isnan;
using std::sqrt;
using std::pow;
#define PI 3.141592653589793

#define bytes2mega (1./1048576.)

#define INTxx int


#include <assert.h>
#include <iostream>
#include "kml/dom.h"
#include "kml/base/file.h"
#include "kml/engine.h"
#include "kml/base/date_time.h"
using kmldom::CoordinatesPtr;
using kmldom::KmlFactory;
using kmldom::PlacemarkPtr;
using kmldom::PointPtr;

using kmldom::FolderPtr;
using kmldom::FeaturePtr;

using kmlbase::Attributes;

string convert_min_to_time(int total_min);

class Footprint3D{
private:
    // ------- Private Functions -----------
	void Destruct_footprint();		// Destructor calls this function
	void Destruct_all_points();		// Destructor calls this function
    
    // Loads the file but not the points
	int load_points_file(string allPointsFileName);
	
	//This is the one that does the actual work, but it must be called from one of the public functions of the same name
	void append_to_existing_footprint();
    
	// --------- Some parameters that should eventually be variables -----------
	const static double arm = 50000;                // Pretty sure this is km
	const static double eps = 1e-14;				//to account for roundoff error when finding 2pi angle 
	
	const static double min_z = 0.;
	const static double max_z = 18.289;	//km
    
    // --------- Parameters ----------------
    // These should be kept and should appear in the copy constructor
	double bin_size;
	int num_bins;
    
    double footprint_delta_t;           //seconds
    int footprint_num_range;
    
    double footprint_launchLat;
    double footprint_launchLon;    
    double footprint_launchAzimuth;     // Degrees  //Azimuth is measured CLOCKWISE from North

    // These are global and may get changed, but never reset.  Actually, delta_t never gets changed.
	double UTC_Initial;
	double UTC_Final;
    
    // Use this to make sure loading at each timestep is used correctly
	int lastTimeStep;   //Check code to make sure this is still useful
    
    // ~~~~~ These should disappear soon ~~~~~~

	// This is number of timesteps over which the points come to us from Architecture class for making a footprint.
	//  This is NOT the number of timesteps associated with any of the propagations
//	int total_num_range;
	
    // I WOULD LIKE TO PHASE OUT THE ALL_POINTS STUFF!!!!
	// For appending points into existing footprints, will need to keep track of some timing information
	double all_points_UTC;
	double all_points_delta_t;          //seconds
	int all_points_num_range;
    
    double all_points_launchLat, all_points_launchLon, all_points_launchAzimuth;    //Azimuth is measured CLOCKWISE from North
	
	double footprint_UTC;       // PHASE THIS OUT
    


	// ---------- Vectors in the order that they appear in the destructor -------------
	vector<vector<vector<vector<Point> > > > footprint_storage3D;	//[cur_index][z_bin][cur_hull][point]
    
    
	vector<vector<Point> > current_points; //[z_bin][point index]
	vector<Point> all_points; //[point index]
	ifstream allPointsFile;
	vector<int> num_all_points_at_tstep;	//This applies to the all_points vector, not to the footprint vector
	vector<vector<Point> > all_points_passed_in;	//[timestep][point index]


	
public:
    // -------------- Accessor Functions -------------------
    double getInitialUTC() { return UTC_Initial; }
//    double getFinalUTC() { return UTC_Final; }
    double getFinalUTC() { return (UTC_Initial + footprint_num_range*footprint_delta_t/(24.*3600.)); }
    double getNumRange() { return footprint_num_range; }
    double getDeltaT() { return footprint_delta_t; }
    double getNumBins() { return num_bins; }
    double getBinSize() { return bin_size; }
    double getLaunchLat() { return footprint_launchLat; }
    double getLaunchLon() { return footprint_launchLon; }
    double getLaunchAzimuth() { return footprint_launchAzimuth; }
    
    void SetFootprintUTC(double curUTC);
    double GetFootprintUTC();
    
    // All in degrees pleeze
    void setAzimuthDeg(double azimuth_in) { footprint_launchAzimuth = azimuth_in; }
    void setLaunchLatDeg(double lat) { footprint_launchLat = lat; }
    void setLaunchLonDeg(double lon) {footprint_launchLon = lon; }

    // ---------- Constructors and Destructors -------------
    Footprint3D();	// Empty constructor
//	Footprint3D(string pointsFile, double bin_size_in, double footprint_UTC_in, double footprint_delta_t_in);
	Footprint3D(string pointsFile, double bin_size_in);
	Footprint3D(string FootprintVectorFile);
    Footprint3D(string SUAFile, double timeOn, double timeOff);
//    Footprint3D(SkyGrid incomingGrid);
    Footprint3D(PointCloud &incoming);
    Footprint3D(PointCloud *incoming, double bin_size_in);
//    Footprint3D(SkyGrid *incoming);


//    Footprint3D(const Footprint3D &incomingFP);
    Footprint3D& operator=( const Footprint3D& rhs );

	~Footprint3D();
    void RemoveFootprintKeepTiming();

	
    // ---------- Functions involved in footprint creation -------------
    // ----- Arranged by order called and indented by heirarchy --------
	void generate_footprint_at_timesteps();
        void destruct_leaves_of_current_points();
        void load_points_at_timestep(int tstep);
        void test_drive_swinging_arm(int timeIX);
            vector<Point> sort_unique_points(int cur_index);
            vector<vector<int> > swing_the_arm(vector <Point> &unique_points);
                int do_they_intersect(vector<Point> *segments);
                void update_available_points(vector<bool> *available, vector<int> *avail_IX, vector<Point> &unique_points,vector<int> &hull_indices);

    // ------------ Functions involved in footprint I/O ----------------
    void store_footprint_as_vector(char *footprintFileName_char);
    void store_footprint_as_vector(string footprintFileName);
	int load_footprint_as_vector(string footprintFileName);
	int load_footprint_as_vector(char *footprintFileName);
	void store_footprint_as_points();                       // Don't think I need this anymore, but keep just in case
    
    // These should get moved to Architecture
	void make_facet_files(string folderName, int startTimeSeconds, int offsetTimeSeconds, int obsolete);    //last argument depreciated and unused
	void make_facet_files(char  *folderName, int startTimeSeconds, int offsetTimeSeconds, int obsolete);
    void exportGoogleEarth(char *googleEarthFile);     // HAS HARDCODED DATE!!!!
    void exportGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min);


    // ----------- Functions involved in footprint merges ---------------
	void append_to_existing_footprint(string newPointsFile);
//	void append_to_existing_footprint(vector<vector<Point> > &total_points_at, double InitialUTC, double tstepMinutes,
//                                      double launchLat, double launchLon, double launchAzimuth);
//	void MergeFootprintVectors(string FP2, string FPOut);
    void MergeFootprintVectors(Footprint3D &incomingFP);
    void MergeFootprintVectors(Footprint3D *incomingFP);    // Cython deals in pointers, not the actual objects themselves
    void SmoothedOut(double newDeltaT = -1);


    // ----------- Functions involved in footprint manipulations ---------------
	void AddToFootprintUTC(double plus_hours, double plus_minutes, string outputFileName = "none");   //Change start time
    // void StretchFootprint()  //Basically clicking the envelope on for awhile
    void ChangeLaunchSiteToDeg(double gdlatIN, double lonIN);
    void ChangeLaunchSiteToRad(double gdlatIN, double lonIN);
	void ChangeAzimuthToDeg(double newAzimuth);
	void ChangeAzimuthToRad(double newAzimuth);

    double ChopTimeAt(double seconds);
    void ShiftFootprintByMinutes(int shiftHowManyMinutes);
    void SlideFootprintBySeconds(int howManySeconds);

    
    vector<vector<Point> > DumpFootprintToPoints();
    vector<Point> CollapseFootprintToPoints();
    int ProjectAllPointsDown();

//    void GridTheSky(vector<vector<Point> > &total_points_at);
//    void ExportBinnedDebrisGoogleEarth();


    
//	void append_to_existing_footprint(string newPointsFile, string existingFootprintFile);

    // GRAVEYARDED
    //	void create_a_footprint(int cur_index);


};



#endif
