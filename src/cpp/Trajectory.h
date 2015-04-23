/*
 *  Trajectory.h
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 5/30/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRAJECTORY_CLASS
#define TRAJECTORY_CLASS

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>

//using namespace boost::numeric::ublas;

//#define ColVec boost::numeric::ublas::vector<double>
//#define MatrixM boost::numeric::ublas::matrix<double>

#include <iomanip>
#include <iostream>
using std::string;
using std::cout;
using std::endl;
using std::cerr;

#include <cmath>
using std::abs;
using std::floor;
using std::pow;
using std::sqrt;
using std::isnan;


#include <fstream>
using std::ofstream;
using std::ifstream;
using std::ios;

#define INTxx int


//// ~~~~~~~~~~~~~~~~ Includes and whatnot for GSL ~~~~~~~~~~~~~~~~~~~~~~~
//#include <gsl/gsl_spline.h>
//#include <gsl/gsl_errno.h>
//
//#include <stdio.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~ END GSL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



// ~~~~~~~~~~~~~~~~ Includes and whatnot for GoogleEarth API ~~~~~~~~~~~~~~~~~~~~~~~
#include <assert.h>     /* assert */
#include "kml/dom.h"  // The KML DOM header.
#include "kml/engine.h"
#include "kml/base/file.h"
#include "kml/base/color32.h"


// libkml types are in the kmldom namespace
using kmldom::CoordinatesPtr;
using kmldom::KmlPtr;
using kmldom::KmlFactory;
using kmldom::PlacemarkPtr;
using kmldom::PointPtr;
using kmldom::MultiGeometryPtr;
// ~~~~~~~~~~~~~~~~~~~~~~~~~ END GoogleEarth API ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//#include <vector>
//using std::vector;

#include "Footprint3D.h"
//#include "Random_Number.h"
//#include "Debris.h"
#include "Point.h"
#include "SkyGrid.h"



#define PI 3.141592653589793

#define TU 806.80415	// sec
#define DU 6378.145		// km  

#define R_equator 6378.1370   //[km]
#define R_polar 6356.7523	  //[km]
#define ecc_Earth 0.081819221456
#define rotEarthRad 0.0000729211585530  //[rad/sec]
#define mu_earth 398600.440	//[km^3/sec^2]




// IMPORTANT NOTES:
//	* Matrix class DOES NOT ALLOW scalar multiplication from the left, only from the right.
//	* Matrix class DOES NOT ALLOW setting two ColVectors equal to each other.  You have to do sth like: v2 = v1 + 0


class Trajectory{
private:
	
//	//Read In
//	double Fmax;
//	double Isp;
//	double Mass0;
//	int NumSteps;
//	
//	double gdlat;
//	double lon;
//	double alt;
//	double Vx0, Vy0, Vz0;
//	
//
	int launch_year;
	int launch_month;
	int launch_day;
	int launch_hours;
	int launch_minutes;
    int launch_seconds;
//	int timezone_offset;
////	double launch_UTC;
//	
//	double min_thrust_time, max_thrust_time;
//	double *ThrustUx;
//	double *ThrustUy;
//	double *ThrustUz;
//	double *ThrustMag;
//	double *ThrustTime;
//	
//	int numAccX, numAccY, numAccZ;
//    double NXstart, NXend, NYstart, NYend, NZstart, NZend;
//
//	int NumStages;
//	double *StagePropMass;
//	double *StageStructMass;
//	double *StageTimes;
//	double *StageISP;
//	double *StageFmax;
//	double PayloadMass;
//	
////	double *RocketStateTimes;
////    int StateTimesLength;
//    
//    vector <double> RocketStateTimes;
//    vector <double> StageDownStateTimes;
//
//
//    int timestepsUntilCutoff;
//	
//	double Initial_UTC;
//	double Aref;
//    double Cd;
//	
//    double launchAzimuth;
//	
//	//Means and Std Devs from GRAM
//	//Assuming UVW is East/North/Up
//	// Altitude is in [km]
//	double min_alt_rho, max_alt_rho, min_alt_wind, max_alt_wind;
//	
//	//Temp is currently coming from someplace else
//	// Temp is in Kelvin, altitude is in km
//	double min_alt_temp, max_alt_temp;
//	
//	
////	// GSL Interpolators
////	gsl_interp_accel *ThrustUx_acc;
////	gsl_spline *ThrustUx_spline;
////	gsl_interp_accel *ThrustUy_acc;
////	gsl_spline *ThrustUy_spline;
////	gsl_interp_accel *ThrustUz_acc;
////	gsl_spline *ThrustUz_spline;
////	gsl_interp_accel *ThrustMag_acc;
////	gsl_spline *ThrustMag_spline;
////	
////	gsl_interp_accel *RhoMean_acc;
////	gsl_spline *RhoMean_spline;
////	gsl_interp_accel *WindUMean_acc;
////	gsl_spline *WindUMean_spline;
////	gsl_interp_accel *WindVMean_acc;
////	gsl_spline *WindVMean_spline;	
////	gsl_interp_accel *WindWMean_acc;
////	gsl_spline *WindWMean_spline;
////	
////	gsl_interp_accel *RhoSd_acc;
////	gsl_spline *RhoSd_spline;
////	gsl_interp_accel *WindUSd_acc;
////	gsl_spline *WindUSd_spline;
////	gsl_interp_accel *WindVSd_acc;
////	gsl_spline *WindVSd_spline;
////	gsl_interp_accel *WindWSd_acc;
////	gsl_spline *WindWSd_spline;
////	
////	gsl_interp_accel *Temp_acc;
////	gsl_spline *Temp_spline;
//	
//	
//	// GSL Gaussian Random Number Generator
////	Random_Number *rand_num;
//	
//	// Uncertain stuff
//	int num_uncert_density_steps;
//	double *Density_uncert_steps;
//	double **Density_profiles;
////	gsl_interp_accel *RhoUncert_acc;
////	gsl_spline *RhoUncert_spline;
//	
//	//Derived
	int StateDot_Option;
//	int Density_Option;
//	int Wind_Option;
	bool Properly_Initialized;
//
//	double time_start, delta_t;
////	int rocket_num_time_steps;
//    
//    int *RocketNumTimeSteps;
//	
//	double initial_radius;
//    double local_elevation;
//
//
//	//These are some variables that will be used within the propagations at each timestep
//	double current_UTC;
//	double current_gdlat;
//	double current_lon;
//	double current_alt;
//	double current_radius_earth;
//	ColVec current_ECEF;
//
//	// These are used for random atmosphere stuff
//	double atm_bin_size_km;
//	int num_wind_bins;
//	std::vector<ColVec> WindProfileUVW;
//
//    double headingAngle0;
//    
//    int debris_state_index;
//
//
//	
////	double last_timestep;
//	static const int state_vec_size = 7;
//	
//	void Get_Nominal_Trajectory(string nominal_file);
//	void Get_Nominal_Suborbital(string nominal_file);
//
//
//
//    
//	void Get_Approximate_Wind(string approx_wind_filename);
//	void Get_Temperature_Data(string temperature_file);
//
//    ColVec DeltaForce;
//	void Get_Random_Thrust_Offsets(int runNumber);
//	ColVec GetThrustForce(double cur_time);
//    
//    double SuborbitalTotalAcc(double cur_time);
//    ColVec GetBodyAcc(double cur_time);
//    ColVec GetBodyAngles(double cur_time);
//    
////	ColVec VwindCurveSEZ(ColVec r__ECI);
//	ColVec VwindCurveUVW(ColVec r__ECI);
//	ColVec VwindCurveUVW();
//
//	double DensityCurve(double altitude);
//	void Initialize_density_profiles(int num_profiles);
//	void load_density_profile(int run_number);
//    
//
//	ColVec CalcVinfECI(double timestep, ColVec State);
//	ColVec CalcAeroDrag(double timestep, ColVec State);
//
//
//	int num_per_batch;
//	std::vector< std::vector< std::vector<double> > > State_Vector_Storage_Vec;	//run cases in batches of 5000
//                                                                    //[batchNum][tstep][state vectors]
//
//    std::vector< std::vector< std::vector<double> > > FallingStageStateVectorStorage;	//run cases in batches of 5000
//                                                                    //[batchNum][tstep][state vectors]
//    
//	// For storing debris before we know how many time steps to keep
//	std::vector< std::vector<double> > State_Vector_Buffer_Vec;	//Single vector, 24*3600 long, each entry is state vector of doubles
//	
//	int write_batch_to_file();  //should get graveyarded
//
//    
//    // ---- Suborbital Stuff
//    vector < double > SuborbitalStateTimes;
//
////    //For Suborbital Pitching Model
////    gsl_interp_accel *PitchAngle_acc;
////	gsl_spline *PitchAngle_spline;
////    gsl_interp_accel *Omega_acc;
////	gsl_spline *Omega_spline;
////    gsl_interp_accel *OmegaDot_acc;
////	gsl_spline *OmegaDot_spline;
//    double timeAngleStart, timeAngleEnd;
//    
//    
//    
//	// Probably don't need both of these as written, a lot of overlap.  TODO: cleanup eventually
//	void write_debris_to_file(unsigned int idNum = 0);
//	void write_debris_to_file(unsigned int idNum, std::vector<std::vector<std::vector<double> > > &LatLonAltStorage);
//
//    // ---- Debris Quantities and Functions
//	double Debris_Mass0;
//	int num_debris_pieces;
////    double *DebrisStateTimes;
////    int debris_max_time_steps;
//
//    vector <double> DebrisStateTimes;
//    
//    int debris_num_per_batch;
//	std::vector< std::vector< std::vector<double> > > Debris_State_Vector_Storage_Vec;	//run cases in batches of 5000
//    
//    std::vector<std::vector<std::vector<double> > > DebrisLatLonAltStorage;
//	
//	// These get read in
//	double DebrisInitialUTC;
//
//    
////    double debris_delta_t;

	
public:
    Trajectory();
    


    // ------> Precomputed lat lon alt stuff
    void loadPrecomputedFile(string inFileName, bool isDegrees, bool isAltMeters);
    void Get_TimeLongLatAlt(string fileName, bool isAltMeters, bool isDegrees);
    
    void loadDebris(vector<double> flatPointArray, vector<int> pointIdArray, vector<int> numTimeSteps, int maxTime, double deltaTsec, double timeOffsetSec, double reactionTimeMinutes);
    
    vector<double> PrecomputedStateTimes;
    std::vector<std::vector<std::vector<double> > > PrecomputedLatLonAltStorage;    //rad, rad, km
    
    void setLaunchTime(int launch_year_in, int launch_month_in, int launch_day_in, int launch_hours_in, int launch_minutes_in, int launch_seconds_in);

    
    int write_to_google_earth_native(string basename, int batch_num, double cutoffAltKM);
    int write_to_google_earth_native(string basename, int printThisMany);
    
//    std::vector<std::vector<std::vector<double> > > TransformToLatLonAlt(int num_to_write);
//    std::vector<std::vector<std::vector<double> > > TransformToLatLonAlt(int num_to_write, double cur_UTC);

    
//    
//    
//    // ---- Things associated with LRHC
//    std::vector< std::vector<Point> > Single_Trajectory_Tube;	//[timestep][point]
//    void SimLeftRightHotCold(int fidelity, double semiMajorPercent, double ecc);
////    void write_LRHC_points(char *outFileName, bool isInstantaneous, double tstepMin, double startUTC, double minutesOn);
//    void SimLeftRightHotCold_Slide();
//
//    
//    void write_LRHC_points_instantaneous(char *outFileName, double tstepMin);
//    void write_LRHC_points_BlockOn(char *outFileName, double tstepMin, double startUTC, double minutesOn);
//    
//    void write_LRHC_points_file(char *outFileName, double tstepMin);
//    void write_stage_down_points_file(char *outFileName, double tstepMin);
//    void write_points_file(char *outFileName, double tstepMin, int StateDot_In);    // This is the function that does all the work
//
//    
//
//    
//    // ---- Things associated with debris propagation / monte carlo
//    void MonteCarloDebris(Footprint3D &my_footprint, string footprintVectorFile);
//    void MonteCarloDebris(Footprint3D &my_footprint, double time_explode_until_in);
//    vector<vector<Point> > assemble_all_points_debris(double tstepMinutes);
//    
//    void MonteCarloDebris(SkyGrid &my_SkyGrid, double time_explode_until_in);
//
//
//    
////	int write_batch_to_google_earth(string basename, int batch_num);
//
//
//    
//    int getTimestepsUntilCutoff();
//
//
//	double get_local_earth_radius(double gdlat /*rad*/);
//
//    void write_trajectory_to_file(char *outputFileName, int num_to_write);
//
////	void write_single_trajectory_to_file(string outputFileName);
////	void write_single_trajectory_to_file(char *outputFileName);
//
//    void read_trajectory_from_file(char *outputFileName);
//
//    
//	double* TestFunction(int s);
//
//
//	
////	Trajectory(char *StateDot_Option_In, char *Wind_Option_In, char *Density_Option_In);
////	Trajectory(string fileName);
//
//    
////	Trajectory();
////	~Trajectory();
//	
//	void Initialize_First_Stage(int num_per_batch_in, double delta_t_in, string nominal_traj_filename);
//	void Initialize_Stages(int num_per_batch_in, double delta_t_in, string nominal_traj_filename);
//	void Initialize_Suborbital(double delta_t_in, string nominal_traj_filename);
//	void Initialize_Single_Debris(double *StateVectorPlusUtcIn, double delta_t_in);
////	std::vector<std::vector<std::vector<double> > > Propagate_Debris_From_Catalog(double *StateVectorPlusUtcIn, Debris CurrentCatalog, double delta_t_in, unsigned int idNum);
////    std::vector<std::vector<std::vector<double> > > Frisco_Propagate_Debris_From_Catalog(double *StateVectorPlusUtcIn, Debris CurrentCatalog, double delta_t_in, unsigned int idNum);
//    
//    
//    
//    std::vector<std::vector<double> > FriscoDebrisPropagation(double *initialState, double* DebrisVel,
//                                                                          double mass, double Sref,
//                                                                          int nCD, double *MinfCD, double *CD,
//                                                                          int CLOption, int nCL, double *MinfCL, double *CL,
//                                                                          double LoverD,
//                                                                          int atmosOption, double *altitudeList, double *densityList, double *UList, double *VList, double *WList,
//                                                                          int GEOptions, char *filename, int nList, int PlanetModel, double dtInterval, int ndtInterval, double thetag0);
//    
//    
//    
//    std::vector<std::vector<double> > Propagate_Stage_Down_To_Ground(double *BreakupStateVector, double BreakupStateTime, int run_number);
//    void Propagate_Stage_Down_To_Ground();
//
//    void PropagateWholeRocket();
//
//	
//	//Some general functions that I'm surprised didn't already exist
//	ColVec Cross_Vectors(ColVec Vec1, ColVec Vec2);
//	double Mod(double x, double y);
//	static double VecNorm(ColVec V);
//	MatrixM Rotation(int axis, double angle);	//Rotate coordinate frame
//
//
//	//Trajectory-specific functions
//	ColVec Geodetic_To_ECEF(double gdlat, double lon, double alt);
//	ColVec ECEF_To_Geodetic(ColVec Vec);
//	ColVec ECEF_To_Geodetic(double x, double y, double z);
//
//	
//	MatrixM ECEF_to_ECI_rotation(double UTC);
//	MatrixM ECI_to_SEZ_rotation(ColVec radius);
//	double Sidereal_Time(double UTC);
//	double UTC_Time(int day, int month, int year, int hours, int minutes, int offset);
//
//	ColVec Statedot_first_stage(double timestep, ColVec State);
//	ColVec Statedot_debris(double timestep, ColVec State);
//	ColVec Statedot_suborbital(double timestep, ColVec State);
//
//
//
////	void Rk4(ColVec state);
//
//	void Propagate_First_Stage();
//    void Propagate_Stages();
//
//	void Propagate_Suborbital();
//	void Propagate_Debris();
//    
//    
//	double HowLongUntilEnterNAS();
//
//	
//	
//	static int func (double t, const double y[], double f[], void *params);
//	static int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
//	int gsl_integrate_first_stage (ColVec State0, int run_number);
//	int gsl_integrate_stages (ColVec State0, int run_number);
//	int gsl_integrate_suborbital (ColVec State0, int run_number);
//	int gsl_integrate_debris (ColVec State0, double *StateTimesBuffer);
//
//	int initialize_random_atmosphere();
//
//	void convert_to_lat_lon();
//
////	int get_num_time_steps();
//
////	void read_debris_from_file();
////	void read_single_trajectory_from_file();
//
//	void Transform_Nominal(string nominalFile, string outputFile);

	
};


#endif
