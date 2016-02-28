/*
 *  Architecture.h
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 7/19/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ARCHITECTURE_CLASS
#define ARCHITECTURE_CLASS

#define PI 3.141592653589793
#define R_equator 6378.1370   //[km]
#define R_polar 6356.7523	  //[km]
#define ecc_Earth 0.081819221456
#include "Trajectory.h"
#include "Debris.h"
#include "Point.h"
#include "Footprint3D.h"

#define INTxx int


#include <cmath>
using std::pow;
using std::sqrt;
using std::sin;
using std::fabs;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ofstream;
using std::ifstream;
using std::ios;

#include <vector>

class Architecture{    
public:
	Architecture();
	
	//These should be merged someday, very similar
	void read_single_trajectory_from_file(char  *trajectoryFileName);
	void read_debris_from_file(unsigned int idNum = 0);
	void reset_debris();

	double get_local_earth_radius(double gdlat /*rad*/);
	void make_tube_around_single_trajectory(double r_desired /*km*/, int fidelity);
	void explode_at_all_points();


//	void write_all_points(string outFileName);
//	void write_all_points_rocket(string outFileName, double tstepMinutes, double debrisCutoff_in);


	//void stitch_together(double tstepMinutes, double debrisCutoff_in);  //GRAVEYARDED
	
//	void make_google_earth_files(string folderName, double debrisCutoff);

//	void SimLeftRightHotCold(char *nominalFileName, int fidelity, double semiMajorPercent, double ecc);
//	void write_LRHC_points(string outFileName, bool isInstantaneous, double tstepMin, double startUTC, double minutesOn);
	vector<vector<Point> > write_all_points_debris(double tstepMinutes, string outFileName = "none");

//	void MonteCarlo();
	void MonteCarlo(char WindOption[], char DensityOption[], string nominal_traj_filename, string footprintVectorFile);

	vector<vector<Point> > ExplodeAtEdges(double tstepMinutes);
	void GoMatlab(vector<vector<Point> > &total_points_at);


private:
	std::vector<std::vector<std::vector<double> > > RocketLatLonAltStorage;
	std::vector<std::vector<std::vector<double> > > DebrisLatLonAltStorage;
	
	// These get read in
	double DebrisInitialUTC;
	int debris_num_time_steps;
	double *DebrisStateTimes;
	
	// This gets calculated
	double DebrisFinalUTC;
	
	double Initial_UTC;
	int num_time_steps;
	double *RocketStateTimes;		// Counts the seconds during launch, t=0 is at Initial_UTC
	int state_vector_size;
	double first_stage_struct_mass;
	
	std::vector< std::vector< std::vector<double> > > State_Vector_Storage_Vec; //[traj][tstep][state_vector_component]
	std::vector< std::vector<Point> > Single_Trajectory_Tube;	//[timestep][point]
	
	
	
};




#endif
