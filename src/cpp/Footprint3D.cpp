/*
 *  Footprint3D.cpp
 *  CompactEnvelope3D
 *
 *  Created by Thomas Colvin on 3/20/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 *  Note: Eventually should add some error checking to this, such as
 *		- In the event of multiple hulls per bin, check that the hulls don't overlap
 *		- When shrinking bin sizes, you can wind up with gaps...check that bin indices are contiguous
 *
 *  * Within these functions, it is appropriate to truncate the points at the top of the NAS.  Nowhere else.
 *
 */

#define PRINT_SUAS false

#define COMPOSITE_MISSION_LAT -365
#define REENTRY_LAT -370

#define R_equator 6378.1370   //[km]
#define R_polar 6356.7523	  //[km]


#include "Footprint3D.h"
//#include "Trajectory.h"


string convert_min_to_time(int total_min);


Footprint3D::Footprint3D() {}

Footprint3D::~Footprint3D(){
//	cout << "In the full destructor" << endl;
	Destruct_footprint();
	Destruct_all_points();
}

void Footprint3D::Destruct_footprint() {
	// footprint_storage3D and boundary_point_storage
	for (int i = 0; i < footprint_num_range; i++) {
		for (int j = 0; j < footprint_storage3D[i].size(); j++) {
			for (int k = 0; k < footprint_storage3D[i][j].size(); k++) {
				footprint_storage3D[i][j][k].clear(); } 
			footprint_storage3D[i][j].clear(); }
		
		footprint_storage3D[i].clear(); }
	footprint_storage3D.clear();
	
	// current_points
	destruct_leaves_of_current_points();
	current_points.clear();
	
	// reset timing info
	footprint_UTC = -1;
	footprint_delta_t = -1;
	footprint_num_range = -1;
	
	// reset altitude bin info
	bin_size = -1;
	num_bins = -1;
    
    // reset launchsite info
    footprint_launchLat = -361;
    footprint_launchLon = -361;
    footprint_launchAzimuth = -361;
	
	return;
}	

void Footprint3D::Destruct_all_points() {
	if (all_points.size() > 0) {
		// all_points (recall this is loaded for each timestep)
		all_points.clear(); }
	
	// num_all_points_at_tstep
	if ( num_all_points_at_tstep.size() > 0 ) { num_all_points_at_tstep.clear(); }
	
	// close the file
	allPointsFile.close();
	
	// the points that were passed in as a vector
	if (all_points_passed_in.size() > 0) {
		all_points_passed_in.clear(); }
	
	// reset lastTimeStep
	lastTimeStep = -1;
	
	// reset timing info
	all_points_UTC = -1;
	all_points_delta_t = -1;
	all_points_num_range = -1;
    
    // reset launchsite info
    all_points_launchLat = -361;
    all_points_launchLon = -361;
    all_points_launchAzimuth = -361;
    
	
	return;
}	

void Footprint3D::destruct_leaves_of_current_points(){
// This is its own function because it gets called often when swinging the arm
	
	if (current_points.size() > 0) {
		for (int i = 0; i < num_bins; i++) {
			if (current_points[i].size() > 0) {
				current_points[i].clear(); } } }
}

Footprint3D::Footprint3D(string FootprintVectorFile) {
	footprint_num_range = -5;	//Set this so it doesn't get accidentally used without being loaded first

//    cout << "\n\n\nIN THE VECTORFILE CONSTRUCTOR ~~~~~~~~~~~~~~~~~~\n\n\n";
	// Load up the file
	load_footprint_as_vector(FootprintVectorFile);

	// Calculate a few things that the main constructor calculates just to be consistent
	UTC_Initial = footprint_UTC;
	UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);

}

void Footprint3D::RemoveFootprintKeepTiming(){
    
    // footprint_storage3D and boundary_point_storage
	for (int i = 0; i < footprint_num_range; i++) {
		for (int j = 0; j < footprint_storage3D[i].size(); j++) {
			for (int k = 0; k < footprint_storage3D[i][j].size(); k++) {
				footprint_storage3D[i][j][k].clear(); }
			footprint_storage3D[i][j].clear(); }
		
//		footprint_storage3D[i].clear();
    }
//	footprint_storage3D.clear();
    
//    footprint_num_range = 0;

    
    return;
}

//Footprint3D::Footprint3D(string pointsFile, double bin_size_in, double footprint_UTC_in, double footprint_delta_t_in) {
//	footprint_UTC = footprint_UTC_in;
//	footprint_delta_t = footprint_delta_t_in;
//	
//	Footprint3D(pointsFile, bin_size_in);
//}

Footprint3D::Footprint3D(string pointsFile, double bin_size_in) {
	// Bin the points by ALTITUDE
	bin_size = bin_size_in;
	num_bins = (INTxx) ceil((max_z - min_z)/bin_size);
	
	//number of bins is fixed, so can allocate this vector of vectors now
	current_points.assign(num_bins,vector<Point>());	
	
	cout << "max_z = " << max_z << " and min_z = " << min_z << endl;
	cout << "numbins = " << num_bins << endl;
	
	// Read in information from the pointsFile and get it ready to be used
	load_points_file(pointsFile);
	
	// Since we haven't generated a footprint yet, the points file will dictate the timing info of the footprint we're about to create
	footprint_num_range = all_points_num_range;
	footprint_UTC = all_points_UTC;
	footprint_delta_t = all_points_delta_t; //seconds
    
    footprint_launchLat = all_points_launchLat;
    footprint_launchLon = all_points_launchLon;
    footprint_launchAzimuth = all_points_launchAzimuth;
	
	UTC_Initial = footprint_UTC;
//	delta_t = footprint_delta_t;
	UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);
	
	// Allocate the footprint vector
	footprint_storage3D.assign(footprint_num_range,vector<vector<vector<Point> > >());
	for (int t = 0; t < footprint_num_range; t++) {
		footprint_storage3D[t].assign(num_bins,vector<vector<Point> >()); }
    
    // Now actually generate the footprint
    generate_footprint_at_timesteps();
    
    // We opened the allPointsFile in load_points_file.  Now it's time to close it.
    allPointsFile.close();
}



Footprint3D::Footprint3D(PointCloud &incomingGrid) {
    
    cout << "\n\n\nIN THE POINTCLOUD CONSTRUCTOR ~~~~~~~~~~~~~~~~~~\n\n\n";

    
    incomingGrid.identifyYourself();
    
    // Bin the points by ALTITUDE
	bin_size = incomingGrid.getZBinHeight();
	num_bins = (INTxx) ceil((max_z - min_z)/bin_size);
    
    all_points_num_range = incomingGrid.getPointsRange();
    all_points_UTC = incomingGrid.getInitialUTC();
    all_points_delta_t = incomingGrid.getDeltaT();  //seconds
    
    all_points_launchLat = incomingGrid.getLaunchLat();
    all_points_launchLon = incomingGrid.getLaunchLon();
    all_points_launchAzimuth = incomingGrid.getLaunchAzimuth();
    
    // Since we haven't generated a footprint yet, the points file will dictate the timing info of the footprint we're about to create
	footprint_num_range = all_points_num_range;
	footprint_UTC = all_points_UTC;
	footprint_delta_t = all_points_delta_t; //seconds
    
    footprint_launchLat = all_points_launchLat;
    footprint_launchLon = all_points_launchLon;
    footprint_launchAzimuth = all_points_launchAzimuth;
    
    UTC_Initial = footprint_UTC;
	UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);

    // Allocate the footprint vector
	footprint_storage3D.assign(footprint_num_range,vector<vector<vector<Point> > >());
	for (int t = 0; t < footprint_num_range; t++) {
		footprint_storage3D[t].assign(num_bins,vector<vector<Point> >()); }
    
    // Snag a copy of all_points_total
    all_points_passed_in = incomingGrid.getAllPoints();
//    cout << "all_points_pass_in size = " << all_points_passed_in.size() << endl;
//    
//    exit(45);
    
    //    num_all_points_at_tstep.assign(all_points_num_range,0);
    //
    //	// read in how many points exist at each time step
    //	int num_here;
    //	for (int t = 0; t < all_points_num_range; t++) {
    //		allPointsFile.read((char *) &num_here, sizeof(num_here));
    //		num_all_points_at_tstep[t] = num_here; }
    
    
//    cout << "allocation done, now generate footprint\n";
    // Now actually generate the footprint
    generate_footprint_at_timesteps();
//    cout << "generated!\n";
//    incomingGrid.generateAllPointsFromKDE(0.4);
    
//    // Properties of the points
//	outfile.write((char *) &num_range, sizeof(num_range));
//	outfile.write((char *) &startUTC, sizeof(startUTC));
//	outfile.write((char *) &tstepMin, sizeof(tstepMin));
//    
//    // Launch site info
//    outfile.write((char *) &gdlat, sizeof(gdlat));
//	outfile.write((char *) &lon, sizeof(lon));
//	outfile.write((char *) &launchAzimuth, sizeof(launchAzimuth));
//    
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//	
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
    
}


// I think Cython can't pass c++ objects, but it can pass pointers
Footprint3D::Footprint3D(PointCloud *incomingGrid, double bin_size_in = -1) {
    
    int objectID = incomingGrid->identifyYourself();
    if (objectID == 1){
        // This is a SkyGrid object
        bin_size = incomingGrid->getZBinHeight();
    }
    else {
        // This is a PointCloud object
        bin_size = bin_size_in; }
    
    if (bin_size < 0){
        cout << "\n\n\nERROR!!!  NEGATIVE BIN SIZE!!!  If you passed in a SKYGRID object, you have to set the z_bin_height!!!!  Expect a segfault.\n\n\n";
        return;
    }
    
    // Bin the points by ALTITUDE
    num_bins = (INTxx) ceil((max_z - min_z)/bin_size);
    
    all_points_num_range = incomingGrid->getPointsRange();
    all_points_UTC = incomingGrid->getInitialUTC();
    all_points_delta_t = incomingGrid->getDeltaT(); //seconds
    
    all_points_launchLat = incomingGrid->getLaunchLat();
    all_points_launchLon = incomingGrid->getLaunchLon();
    all_points_launchAzimuth = incomingGrid->getLaunchAzimuth();
    
    // Since we haven't generated a footprint yet, the points file will dictate the timing info of the footprint we're about to create
	footprint_num_range = all_points_num_range;
	footprint_UTC = all_points_UTC;
	footprint_delta_t = all_points_delta_t; //seconds
    
    footprint_launchLat = all_points_launchLat;
    footprint_launchLon = all_points_launchLon;
    footprint_launchAzimuth = all_points_launchAzimuth;
    
    UTC_Initial = footprint_UTC;
	UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);

    // Allocate the footprint vector
	footprint_storage3D.assign(footprint_num_range,vector<vector<vector<Point> > >());

	for (int t = 0; t < footprint_num_range; t++) {
		footprint_storage3D[t].assign(num_bins,vector<vector<Point> >());
    }
        
    // Snag a copy of all_points_total
    all_points_passed_in = incomingGrid->getAllPoints();
    
    // Now actually generate the footprint
    generate_footprint_at_timesteps();


    
    // DEBUG?
    bool debugHere = false;    
    if (debugHere){
        int footprint_tstep = 15;
        int zbinIX = 0;
        int num_hulls_here = footprint_storage3D[footprint_tstep][zbinIX].size();
        cout << "%footprint_tstep = " << footprint_tstep << endl;
        
        for (int hullIX = 0; hullIX < num_hulls_here; hullIX++) {
            cout << "hull{" << hullIX+1 << "} = [\n";
            int num_points_here = (INTxx) footprint_storage3D[footprint_tstep][zbinIX][hullIX].size();
            
            for (int pt = 0; pt < num_points_here; pt++) {
                cout << std::setprecision(12);
                cout << "  " << footprint_storage3D[footprint_tstep][zbinIX][hullIX][pt].get_x()
                << "  " << footprint_storage3D[footprint_tstep][zbinIX][hullIX][pt].get_y() << endl;
            }
            cout << "];\n\n";
        }
        
        cout << "plot(hull{1}(:,1), hull{1}(:,2))\n\n\n";
    }
    
}


Footprint3D& Footprint3D::operator=( const Footprint3D& rhs ) {
    bin_size = rhs.bin_size;
	num_bins = rhs.num_bins;
    
	footprint_num_range = rhs.footprint_num_range;
    footprint_UTC = rhs.footprint_UTC;
	footprint_delta_t = rhs.footprint_delta_t;           //seconds
    
    footprint_launchLat = rhs.footprint_launchLat;
    footprint_launchLon = rhs.footprint_launchLon;
    footprint_launchAzimuth = rhs.footprint_launchAzimuth;


    UTC_Initial = rhs.UTC_Initial;
	UTC_Final = rhs.UTC_Final;
    
    //launchLon
    
    footprint_storage3D = rhs.footprint_storage3D;
    
    return *this;
}






/*! This function does not actually load the points into all_points, but rather it opens the file and reads some of the necessary
 information about the contents.  The all_points vector is populated at each time step by load_points_at_timestep()
 Currently can only load one file at a time; no appending multiple files together

 ---- ALL Points files should have this format -----
 int			number_of_time_steps
 double			UTC_at_start
 double			delta_t_in_minutes
 ints			for (number_of_time_steps) { number_of_points_at_this_time_step }
 vector<Point>	for (number_of_time_steps) { vector_of_points_at_this_time_step }     
 */

int Footprint3D::load_points_file(string allPointsFileName) {
    // Make sure we're dealing with a clean slate
	Destruct_all_points();
	
	// Quick error check
	if (allPointsFile.is_open()) {
		cout << "ERROR: allPointsFile already open.  Please close the old points file manually to be sure this is what you want to do~~~~~~" << endl; 
		//exit(232734);	
	}
	
	// Open the file
	allPointsFile.open(allPointsFileName.c_str(), ios::in | ios::binary);
	if (!allPointsFile.is_open()) {
		cout << "ERROR: failed to open " << allPointsFileName.c_str() << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		exit(53628920);	}
		
	// Read in timing information
	allPointsFile.read((char *) &all_points_num_range, sizeof(all_points_num_range));	
	allPointsFile.read((char *) &all_points_UTC, sizeof(all_points_UTC));
	allPointsFile.read((char *) &all_points_delta_t, sizeof(all_points_delta_t));
    
    // Read in launch site information
    allPointsFile.read((char *) &all_points_launchLat, sizeof(all_points_launchLat));
	allPointsFile.read((char *) &all_points_launchLon, sizeof(all_points_launchLon));
	allPointsFile.read((char *) &all_points_launchAzimuth, sizeof(all_points_launchAzimuth));

    cout << "DEBUG" << endl;
    cout << all_points_launchLat << "   " << all_points_num_range << "   " << endl;
    cout << all_points_launchLon << "   " << all_points_UTC << "   " << endl;
    cout << all_points_launchAzimuth << "   " << all_points_delta_t << "   " << endl;
		
	num_all_points_at_tstep.assign(all_points_num_range,0);
	
	// read in how many points exist at each time step
	int num_here;
	for (int t = 0; t < all_points_num_range; t++) {
		allPointsFile.read((char *) &num_here, sizeof(num_here));	
		num_all_points_at_tstep[t] = num_here; }
	
	
	// Used for debugging purposes only
	bool debugnow = false;
	if (debugnow) {
		static int loadCounter = 1;
		cout << "Loading Points: " << loadCounter << endl;
		loadCounter++;
		
		string debugFileName("GeneratedFiles/LOADPOINTS_PLAINTEXT_DEBUG.txt");
		ofstream debugfile;
		debugfile.open(debugFileName.c_str(), ios::out);
		
		debugfile << all_points_num_range << endl << all_points_UTC << endl << all_points_delta_t << endl;
		
		for (int t = 0; t < all_points_num_range; t++) {
			debugfile << "num_all_points_at_tstep[" << t << "] = " << num_all_points_at_tstep[t] << endl; }
		
		debugfile.close();
	}
	// End of debugging code
	
	
	return all_points_num_range;
}


void Footprint3D::load_points_at_timestep(int tstep) {
    // There are two ways to call this function.  Either you have:
    //   1.) Already set up a file to be read in with the function load_points_file
    //   or
    //   2.) You have loaded some points into all_points_passed_in via the function append_to_existing_footprint
 
	// If we read in a file of points, then the allPointsFile will be open
	if (allPointsFile.is_open()) {
		if (lastTimeStep +1 == tstep) {
			all_points.clear();
			all_points.assign(num_all_points_at_tstep[tstep],Point());
			
			allPointsFile.read((char *) &all_points[0], all_points.size()*sizeof(Point));
			lastTimeStep = tstep;
		}
		else {
			cout << "ERROR!!!! YOURE ABUSING THE load_points_at_timestep FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			exit(112930);
		}
	}
	// If the file is not open, then that means we were directly passed a vector containing the points
	else {
		all_points.clear();
		all_points = all_points_passed_in[tstep];
	}

	return;
}


// Already have a bunch of points 
void Footprint3D::generate_footprint_at_timesteps() {
    
	//pretty sure this should be footprint_num_range and not all_points_num_range
    for (int i = 0; i < footprint_num_range; i++) {				
        // Clean things up
//        cout << "destruct\n";
		destruct_leaves_of_current_points();
		
//        cout << "load" << endl;
		// Load up the points from the current timestep into all_points
		load_points_at_timestep(i);
				
//        cout << "drive" << endl;
        // Make the footprint at this timestep
		test_drive_swinging_arm(i);
	}
}




// ---- ALL Footprint Vector files should have this format -----
// int				number_of_time_steps
// double			UTC_at_start
// double			delta_t_in_minutes
// int				size_of_bins_in_km		NOTE: This is constant, no need to loop on it later
// ints				Loop through the sizes of the vector structure
// vector<Point>	Loop through the vector_of_points at_this_timestep, in_this_altitude_bin, in_this_hull
//
// Note: Vector structure is [time_step][z_bin][cur_hull][point]

void Footprint3D::store_footprint_as_vector(char *footprintFileName_char) {
    string footprintFileName(footprintFileName_char);
    return store_footprint_as_vector(footprintFileName);
}


void Footprint3D::store_footprint_as_vector(string footprintFileName) {
	// In addition to the vector itself, will need some extra information like delta_t, delta_z, and arm length
	//   At the moment, these things are hardcoded
	//must also save the timing information found in points files
	
	ofstream outfile;
	outfile.open(footprintFileName.c_str(), ios::out | ios::binary);

	// Write the timing information
	outfile.write((char *) &footprint_num_range, sizeof(footprint_num_range));
	outfile.write((char *) &footprint_UTC, sizeof(footprint_UTC));
	outfile.write((char *) &footprint_delta_t, sizeof(footprint_delta_t));
	
	// Write the altitude bin information
	outfile.write((char *) &bin_size, sizeof(bin_size));
	outfile.write((char *) &num_bins, sizeof(num_bins));
    
    // Write launch site information
    outfile.write((char *) &footprint_launchLat, sizeof(footprint_launchLat));
	outfile.write((char *) &footprint_launchLon, sizeof(footprint_launchLon));
    outfile.write((char *) &footprint_launchAzimuth, sizeof(footprint_launchAzimuth));

	
	// Write the sizes of the the vectors and subvectors
	for (int t = 0; t < footprint_num_range; t++) {
		for (int ix = 0; ix < num_bins; ix++) {
			int num_hulls = (INTxx) footprint_storage3D[t][ix].size();
			outfile.write((char *) &num_hulls, sizeof(num_hulls));

			int num_points;
			for (int jx = 0; jx < num_hulls; jx++) {
				num_points = (INTxx) footprint_storage3D[t][ix][jx].size();
				outfile.write((char *) &num_points, sizeof(num_points)); }

			// What if there were no hulls?  Do this just so that we're never NOT writing something
			//  I don't think this is strictly necessary, and it's a little ugly, but it makes me feel better
			if (num_hulls == 0) {
				num_points = 0;
				outfile.write((char *) &num_points, sizeof(num_points)); }
		}}

	// Now write the data
	Point EmptyPoint;
	for (int t = 0; t < footprint_num_range; t++) {
		for (int ix = 0; ix < num_bins; ix++) {
			int num_hulls = (INTxx) footprint_storage3D[t][ix].size();

			for (int jx = 0; jx < num_hulls; jx++) {
				int num_points = (INTxx) footprint_storage3D[t][ix][jx].size();
				outfile.write((char *) &footprint_storage3D[t][ix][jx][0], num_points*sizeof(Point));
				for (int kx = 0; kx < num_points; kx++) {
					//cout << footprint_storage3D[t][ix][jx][kx] << endl;
				}
			}

			//Again, trying to avoid the situation where we don't write anything
			//  So write a single point with the understanding that we are going to throw it out when reading in
			if (num_hulls == 0) {
				outfile.write((char *) &EmptyPoint, sizeof(Point)); }

		}}

	outfile.close();
	
	
	return;
}

int Footprint3D::load_footprint_as_vector(char *footprintFileName_char) {
    string footprintFileName(footprintFileName_char);
    return load_footprint_as_vector(footprintFileName);
}

int Footprint3D::load_footprint_as_vector(string footprintFileName) {
//NOTE: BE CAREFUL USING THIS FUNCTION!!!  There is no current plan for updating the universal timing info i.e. UTC_Initial
//  based on this, so if you try using that variable after loading this up...well...performance not guaranteed.
//  It is likely that the global timing info shoudl just be set here in this function to be equal to the footprint timing info
	
	// Things that get read in
	//  * footprint_num_range
	//  * footprint_UTC
	//  * footprint_delta_t
	//  * bin_size
	//  * num_bins
	//  * footprint_storage3D
    

    

	Destruct_footprint();
	
	ifstream inFile;
	inFile.open(footprintFileName.c_str(), ios::in | ios::binary);
	

	if (inFile.fail()) {
		cout << "Failed to open " << footprintFileName << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		exit(9347590); }

	// Read timing information
	inFile.read((char *) &footprint_num_range, sizeof(footprint_num_range));
	inFile.read((char *) &footprint_UTC, sizeof(footprint_UTC));
	inFile.read((char *) &footprint_delta_t, sizeof(footprint_delta_t));
	
	// Read altitude bin information
	inFile.read((char *) &bin_size, sizeof(bin_size));
	inFile.read((char *) &num_bins, sizeof(num_bins));
    
    // Read launch site information
    inFile.read((char *) &footprint_launchLat, sizeof(footprint_launchLat));
	inFile.read((char *) &footprint_launchLon, sizeof(footprint_launchLon));
	inFile.read((char *) &footprint_launchAzimuth, sizeof(footprint_launchAzimuth));

	// Allocate / read in the structure and information of the footprint vector
	footprint_storage3D.assign(footprint_num_range, vector<vector<vector<Point> > >());
	
	for (int t = 0; t < footprint_num_range; t++) {
		footprint_storage3D[t].assign(num_bins, vector<vector<Point> >());

		for (int ix = 0; ix < num_bins; ix++) {
			int num_hulls;
			inFile.read((char *) &num_hulls, sizeof(num_hulls));
			footprint_storage3D[t][ix].assign(num_hulls, vector<Point>());

			int num_points;
			for (int jx = 0; jx < num_hulls; jx++) {
				inFile.read((char *) &num_points, sizeof(num_points));
				footprint_storage3D[t][ix][jx].assign(num_points, Point()); }

			//Gotta read in that trivial num_points and throw it away
			if (num_hulls == 0) {
				inFile.read((char *) &num_points, sizeof(num_points)); }
		}}
	
	Point GarbagePoint;

	// Now read the data itself
	for (int t = 0; t < footprint_num_range; t++) {
		for (int ix = 0; ix < num_bins; ix++) {
			int num_hulls = (INTxx) footprint_storage3D[t][ix].size();

			for (int jx = 0; jx < num_hulls; jx++) {
				int num_points = (INTxx) footprint_storage3D[t][ix][jx].size();
				inFile.read((char *) &footprint_storage3D[t][ix][jx][0], num_points*sizeof(Point));
//                cout << "this Point " << footprint_storage3D[t][ix][jx][0] << endl;

            }
            

			//Gotta read in that empty point and throw it away
			if (num_hulls == 0) {
				inFile.read((char *) &GarbagePoint, sizeof(GarbagePoint));
			}
		}}


	inFile.close();
	
	return footprint_num_range;
}


void Footprint3D::make_facet_files(char *folderName, int startTimeSeconds, int offsetTimeSeconds, int obsolete){
    string folderNameStr(folderName);
    make_facet_files( folderNameStr,  startTimeSeconds,  offsetTimeSeconds,  obsolete);
    return;
}

void Footprint3D::make_facet_files(string folderName, int startTimeSeconds, int offsetTimeSeconds, int obsolete){
	int base_time = startTimeSeconds - offsetTimeSeconds;	// The time when the space mission begins.  Might not have an envelope if debris is high enough.
	int delta_time = footprint_delta_t;

	// Use this one for generating stuff that will run on the iMac at Ames
    //	string pathToPreferenceFiles = "/Users/tcolvin1/Documents/workspace3p6/facet/user/preferences/";
    //	string pathToPreferenceFiles = "ToCopy/";
    
    // We're going to create a folder (myFolder) and place it in the "work" directory of the FACET code tree.
    // Facet/user/work/ScenarioName/VehicleMissionName/ will contain all the files needed to run a specific instance
    //  If we're consistent about this, then we can always say that the preference files are at 
    // ../../../preferences/
    
    string pathToPreferenceFiles = "../../../preferences/";


	if (offsetTimeSeconds < 0){
		cout << "You cannot have a negative offset time!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		exit(23948); }
	
	string suaFileName("SUA_GeneratedSUA");
	
//	cout << "folderName = " << folderName + suaFileName << endl;

	// GENERATE AN SUA FILE FOR FACET
	// Need to convert into lat/lon as well as change units to ft
	double KmToFt = 3280.8399;
	ofstream SUA;
	SUA.open( (suaFileName).c_str() );
//	SUA.open( (folderName + suaFileName).c_str() );
//    if (SUA.is_open()){
//		cout << folderName + suaFileName + "  SUA file opened successfully!\n";
//	}
//	else {
//		std::cerr << folderName + suaFileName + "  SUA file failed to open :(\n";	}
    
    
	double alt_base = 0;
	int counter = 0;
	// t = timestep, z = z_bin, h = hull, i = a point
	for (int t = 0; t < footprint_num_range; t++) {
		int start_time = base_time + t*delta_time;
		int stop_time = base_time + (t+1)*delta_time;
		
		int curr_num_bins = (INTxx) footprint_storage3D[t].size();
		for (int z=0; z < curr_num_bins; z++){
			int cur_alt_lo = (INTxx) floor((alt_base + bin_size*z)*KmToFt);
			int cur_alt_hi = (INTxx) floor((alt_base + bin_size*(z+1))*KmToFt);
			
			int curr_num_hulls = (INTxx) footprint_storage3D[t][z].size();
			for (int h=0; h < curr_num_hulls; h++){
				SUA << "SUA" << t << "_" << z << "\n\n0\n";
				SUA << cur_alt_lo << " " << cur_alt_hi << " " << start_time << " " << stop_time << endl;	//wants time in seconds
				SUA << "0 0 0 1" << endl;
				SUA << footprint_storage3D[t][z][h].size() << endl;
				
				int curr_num_points = (INTxx) footprint_storage3D[t][z][h].size();
				for (int i = 0; i < curr_num_points; i++) {
//					double xm = footprint_storage3D[t][z][h][i].get_x();
//					double ym = footprint_storage3D[t][z][h][i].get_y();
//					
//					//when using oblateness, will need to store the local_R in order to find the correct gdlat values here
//					double Lat = (180*ym)/(PI*footprint_storage3D[t][z][h][i].get_R_local());
//					double Lon = (180*xm)/(PI*R_equator);
                    
                    double Lat = footprint_storage3D[t][z][h][i].get_gdLatDeg();
                    double Lon = footprint_storage3D[t][z][h][i].get_lonDeg();
                    
					
					if (Lon < 0) {
						Lon = 360 + Lon; }
					
					SUA << Lat << " " << Lon << endl; 
				} 			
				counter++;
			}
			
		}
	}
	
	SUA.close();
	
	
	
	string filter_label ("SUA_GeneratedSUA");
	string filter_out_file = "OutfileAcCount";
	
	// GENERATE FILTER FILES FOR FACET THAT USE THE SUAs JUST CREATED
	
	//http://stackoverflow.com/questions/303371/whats-the-easiest-way-to-generate-xml-in-c
	// Has list of ranked xml libraries for c++ in case it ever seems like i need to go that route
	
	string acFilterPrefs("GeneratedSUA_aircraft_filter_preferences.xml");
	string acFilter("GeneratedSUA_aircraft_filters.xml");
	string plotDataPrefs("GeneratedSUA_plot_data_preferences.xml");
	string mainPrefs("GeneratedSUA_MainPreferences.xml");
	string acPrefs("GeneratedSUA_aircraft_preferences.xml");
	string airspace("GeneratedSUA_airspace.xml");
	string airspacePrefs("GeneratedSUA_airspace_preferences.xml");
	string cwamPrefs("GeneratedSUA_cwam_preferences.xml");
	string toCopyFolder("ToCopy/");
	
	string cp("cp ");
	string mv("mv ");

	
	
	ofstream filter;
	filter.open((acFilterPrefs).c_str());
	
	filter << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl << endl;
	filter << "<aircraft_filter_preferences>" << endl;
	
	filter << "\t<aircraft_filter>" << endl;
	filter << "\t\t<fca_filter>";	
	for (int t = 0; t < footprint_num_range; t++) {
		for (unsigned int z = 0; z < footprint_storage3D[t].size(); z++) {
			filter << "SUA" << t << "_" << z << ", "; } }
	
	//We wrote one comma and space too many, must back up a two characters and overwrite them
	long pos;
	pos = filter.tellp();
	filter.seekp(pos-2);
	filter << "</fca_filter>" << endl;

	
	filter << "\t\t<user_filter_label>";	
	for (int t = 0; t < footprint_num_range; t++) {
		for (unsigned int z = 0; z < footprint_storage3D[t].size(); z++) {
			filter << "SUA" << t << "_" << z << ", "; } }
	
	//Wrote too many, must back up a two characters and overwrite
	pos = filter.tellp();
	filter.seekp(pos-2);
	filter << "</user_filter_label>" << endl;
	
	filter << "\t\t<filter_color>0,255,255</filter_color>" << endl;
	filter << "\t\t<aircraft_symbol>dots</aircraft_symbol>" << endl;
	filter << "\t\t<show_aircraft/>" << endl;
	filter << "\t\t<show_histories/>" << endl;
	filter << "\t</aircraft_filter>" << endl;	
	
	filter << "</aircraft_filter_preferences>" << endl;
	filter.close();
	
	
	
	// Make the other filter file
	filter.open( (acFilter).c_str() );
	filter << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	filter << "<!DOCTYPE aircraft_filters SYSTEM \"" + pathToPreferenceFiles + "aircraft_filters.dtd\" [" << endl;
//	filter << "<!DOCTYPE aircraft_filters SYSTEM \"../../../../preferences/aircraft_filters.dtd\" [" << endl;
	filter << "\t<!ENTITY aircraft_filters SYSTEM \"GeneratedSUA_aircraft_filter_preferences.xml\">" << endl;
	filter << "]>" << endl;
	
	filter << "<aircraft_filters>" << endl;
	filter << "\t&aircraft_filters;" << endl;
	filter << "</aircraft_filters>" << endl;
	filter.close();
	
	
	
	
	filter.open( (plotDataPrefs).c_str() );
	filter << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl << endl;
	
	filter << "<plot_data_preferences>" << endl;
	filter << "\t<plot_window>" << endl;
	filter << "\t\t<plot_calculations>" << endl;
	filter << "\t\t\t<plot_number_of_ac/>" << endl;
	filter << "\t\t\t<for_filtered_ac>" << endl;
	
	filter << "\t\t<filter_name>";	
	for (int t = 0; t < footprint_num_range; t++) {
		for (int z = 0; z < footprint_storage3D[t].size(); z++) {
			filter << "SUA" << t << "_" << z << ", "; } }
	
	//Wrote too many, must back up a two characters and overwrite
	pos = filter.tellp();
	filter.seekp(pos-2);
	filter << "</filter_name>" << endl;
	
	filter << "\t\t\t</for_filtered_ac>" << endl;
	filter << "\t\t\t<where_flying/>" << endl;
	filter << "\t\t\t<output_filename>" << filter_out_file << "</output_filename>" << endl;
	filter << "\t\t</plot_calculations>" << endl;
	filter << "\t</plot_window>" << endl;
	

	filter << "</plot_data_preferences>" << endl;
	filter.close();
	
	
	filter.open( (mainPrefs).c_str() );
	filter << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl;
	filter << "<!DOCTYPE preferences SYSTEM \"" + pathToPreferenceFiles + "preferences.dtd\" [" << endl;

//	filter << "<!DOCTYPE preferences SYSTEM \"../../../../preferences/preferences.dtd\" [" << endl;
	filter << "<!ENTITY aircraft_filter_preferences SYSTEM \"GeneratedSUA_aircraft_filter_preferences.xml\">" << endl;
	filter << "<!ENTITY aircraft_preferences SYSTEM \"GeneratedSUA_aircraft_preferences.xml\">" << endl;
	filter << "<!ENTITY airspace_preferences SYSTEM \"GeneratedSUA_airspace_preferences.xml\">" << endl;	
	filter << "<!ENTITY plot_data_preferences SYSTEM \"GeneratedSUA_plot_data_preferences.xml\">" << endl;
	filter << "<!ENTITY cwam_preferences SYSTEM \"GeneratedSUA_cwam_preferences.xml\">" << endl;			//don't have these at the moment
	filter << "]>" << endl;
	
	filter << "<preferences>" << endl;
	filter << "&aircraft_preferences;" << endl;
	filter << "&airspace_preferences;" << endl;
	filter << "&plot_data_preferences;" << endl;
	filter << "&cwam_preferences;" << endl;
	filter << "</preferences>" << endl;
	filter.close();
	
	
	filter.open( (acPrefs).c_str() );
	filter << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	filter << "<aircraft_preferences>" << endl;
	filter << "<!-- triangles, images, dots, dots_with_heading, pixels -->" << endl;
	filter << "<display_symbols>pixels</display_symbols>" << endl;
	filter << "<display_identical_symbols/>" << endl;
	filter << "<aircraft_tag>" << endl;
	filter << "<display_id_tag/>" << endl;
	filter << "</aircraft_tag>" << endl;
	filter << "<show_all_aircraft/>" << endl;
	filter << "<show_mouse_over_data_blocks/>" << endl;
	filter << "<!-- Aircraft filters -->" << endl;
	filter << "&aircraft_filter_preferences;" << endl;
	filter << "</aircraft_preferences>" << endl;
	filter.close();
	
	filter.open( (airspace).c_str() );
	filter << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	filter << "<!DOCTYPE airspace SYSTEM \"" + pathToPreferenceFiles + "airspace.dtd\" [" << endl;

//	filter << "<!DOCTYPE airspace SYSTEM \"../../../../preferences/airspace.dtd\" [" << endl;
	filter << "<!ENTITY airspace SYSTEM \"GeneratedSUA_airspace_preferences.xml\">" << endl;
	filter << "]>" << endl << endl;
	
	filter << "<airspace>" << endl;
	filter << "&airspace;" << endl;
	filter << "</airspace>" << endl;
	filter.close();
	
	
	//NOTE: Probably not worth actually making the airspace_preferences file.  Just generate one in FACET and use that.
	//  Should probably use system() to copy the file into the current folder 
	
	// Copy the files into the appropriate directory
    system( ("mkdir " + folderName).c_str() );
	system( (mv + suaFileName + " " + folderName).c_str());
	system( (mv + acFilterPrefs + " " + folderName).c_str() );
	system( (mv + acFilter + " " + folderName).c_str() );
	system( (mv + plotDataPrefs + " " + folderName).c_str() );
	system( (mv + mainPrefs + " " + folderName).c_str() );
	system( (mv + acPrefs + " " + folderName).c_str() );
	system( (mv + airspace + " " + folderName).c_str() );
	system( (cp + toCopyFolder + airspacePrefs + " " + folderName).c_str() );
	system( (cp + toCopyFolder + cwamPrefs + " " + folderName).c_str() );
	

	return;
}


// TODO: No longer using minutes, but this could be useful in seconds.
string convert_min_to_time(int total_min) {
	int hours = (INTxx) floor(((double) total_min)/60.0);
	int min = total_min%60;
	
	char buffer[5];
	sprintf(buffer, "%.2d:%.2d", hours,min);
	
	string ans(buffer);
	return ans;
	
}






vector<Point> Footprint3D::sort_unique_points(int cur_bin){
    
    
    // Start by making a list (vec should work too?) out of the current points
    vector<Point> points = current_points[cur_bin];
    int numCurPts = points.size();
    
    // Round the points to the nearest 0.001(?) and put them in a map of maps (which will order them)    
    map<double, map<double, double> > mapToSort;
    map<double, map<double, double> >::reverse_iterator ity;
    map<double,double>::reverse_iterator itx;
    double trunc = 1e-8;
    
    std::pair<map<double,double>::iterator,bool> retVal;
    int num_remaining_pts = 0;
    
    for (int pt = 0; pt < numCurPts; pt++){
        double x = points[pt].get_x();
        double y = points[pt].get_y();
        
        // Since we're already in a zbin, we don't even care what the z-value is anymore.  Could eventually go down to a single map?
        //  Duplicates will be automatically sorted out based on the truncation length.
//        mapToSort[y-fmod(y,trunc)][x - fmod(x,trunc)] = points[pt].get_z();
        
        // Wnat to know how many unique points there are...
        //  This should HOPEFULLY
        //      1.) Create or Access the y-value of the map
        //      2.) Attempt to insert the pair of values for the x and z
        //      3.) If insertion fails because the x-val already exists, then retVal.second = false
        retVal = mapToSort[y-fmod(y,trunc)].insert(std::pair<double,double> (x - fmod(x,trunc), points[pt].get_z()));
        
        if (retVal.second == true){
            // This was a unique point
            ++num_remaining_pts;
        }
    
    }
    
//    cout << "numCurPts = " << numCurPts << "   num_remaining = " << num_remaining_pts <<  endl;
    
    Point tempPt;
    
    vector<Point> remaining_points;
	remaining_points.assign(num_remaining_pts,Point());

//    cout << "ptsHere = [\n";
    int curIX = 0;
    for (ity = mapToSort.rbegin(); ity != mapToSort.rend(); ++ity){
        double yindex = ity->first;

        for (itx = mapToSort[yindex].rbegin(); itx != mapToSort[yindex].rend(); ++itx){
            double xindex = itx->first;
            
            tempPt.set_xyz(xindex, yindex, itx->second);
            remaining_points[curIX] = tempPt;
            
//            cout << xindex << "   " << yindex << endl;
            
//            cout << "remain = " << remaining_points[curIX] << endl;
            
//            if (curIX != 0){
//                double dist = remaining_points[curIX].calc_dist_xy(remaining_points[curIX-1]);
////                if (dist < 1e-1){
//                    cout << "NOOOOOOO " << dist << endl;
////                }
//            }
            
//            cout << "  " << yindex << "  " << xindex << endl;
            ++curIX;
        }
    }
//    cout << "];\n\n\n\n";
    
//	int num_remaining_pts = (INTxx) mylist.size();
    
    return remaining_points;

    
    
    // Look for duplicates
    
    // Load the map into the return vector
    
//    // Start by making a list (vec should work too?) out of the current points
//    vector<Point> points = current_points[cur_bin];
//    
//    // Set the sortBy flag in the point comparator so that we'll get proper behavior
//    
//    // Run the sort
//    
//    // Loop through sorted list, comparing points to neighbors, removing duplicates
//    
//    // Load the list into the return vector (may need to swap / iterate backwards)
    
    
    
    
    
    
	
//	//Find number of points in current bin and allocate space to store them in a working vector 'points'
//	int vec_size = (INTxx) current_points[cur_bin].size();
//	vector<Point> points;
//	points.assign(vec_size,Point());			//allocate the space needed as origin pts.  ALL will get overwritten.
//	
//	//If there are no points, return an empty vector
//	if (vec_size == 0) {
//		return points;	}
//	
//	//Load up the points that are in the current bin
//	for (int i = 0; i < vec_size; i++) {
//		points[i] = current_points[cur_bin][i];
//        //		if (points[i].get_y() < 1e-6) {
//        //			cout << "Shit... = " << points[i].get_y() << endl; }
//		points[i].set_id(i);
//        //		cout << "point" << i << " = " << points[i] << endl;
//        
//	}
//	
//	//sort points by y value with a list
//	list<Point> mylist;
//	list<Point>::iterator it;
//	
//	bool need_origin = true;		//What's going on here?
//	int list_size;
//	it = mylist.begin();
//	
//	int first_index = 0;
//	
//	//Automatically accept the first point
//	mylist.insert(it, points[first_index]);
//	
//	for (int i = first_index+1; i < vec_size; i++) {
//		
//		bool duplicate = false;
//		
//		list_size = (INTxx) mylist.size();
//		int j=0;	//makes sure we don't overflow the size of the present list
//		
//		double current_y = points[i].get_y();
//		double current_x = points[i].get_x();
//		
//		//check to see if the current pt is effectively a duplicate of another point already seen
//		it = mylist.begin();
//		for (int ix = 0; ix < list_size; ix++) {
//			double dist_between = sqrt(pow(current_y - it->get_y(),2) + pow(current_x - it->get_x(),2));
//			if (dist_between < 0.001) {
//				duplicate = true;
//				break;
//			}
//			it++; }
//		
//		if (!duplicate) {
//			//reset the iterator and find where the current point should fit in the list
//			it = mylist.begin();
//			while ((current_y <= it->get_y()) && (j < list_size)) { //finding where to place new point
//				
//				if (current_y == it->get_y()) {
//					//				cout << "TWO POINTS SHARE SAME YVAL!!!!!" << endl;
//				}
//				
//				j++;
//				it++;
//			}
//		}
//		
//		
//		//BIG KLUDGE!!!! FIX THIS SOMETIME!!!
//		//having y=0 x>0 multiple points messes up algorithm, need to sort pts by xval for degenerate yvals
//		//until then, just zero out the x-parts for points with y=0
//		//overload comparison operators in point class.
//		
//		//This is still not totally right...sorting can still go awry near the end of the sorted list...
//		if ((points[i].get_y() == 0) && need_origin && !duplicate) {
//            cout << "I dont think this will ever realistically happen.  ERROR IF IT DOES THOUGH!!!!!!!!!~~~~~~~~~~~\n\n\n\n\n\n\n\n\n";
//            cout << points[i].get_y() << " " << need_origin << " " << duplicate << endl;
//            //			points[i].set_x(0);
//            //			mylist.insert(it, points[i]);  //insert it into the list
//            //			need_origin = false;
//        }
//		else if ((list_size != 1) && (j+1 == list_size) && !duplicate) {	//put point at end of list
//			mylist.push_back(points[i]); }
//		else if (!duplicate) {				//Why do i care if it's positive????~~~~~~~~~~~~
//			mylist.insert(it, points[i]);  //insert it into the list (before element at 'it')
//		}
//		else {}
//	}
//	
//	it = mylist.begin();
//	vector<Point> remaining_points;
//	int num_remaining_pts = (INTxx) mylist.size();
//	remaining_points.assign(num_remaining_pts,Point());
//	
//	for (int i = 0; i < num_remaining_pts; i++) {
//		remaining_points[i] = *it;
//		it++; }
//	return remaining_points;
	
}





//// ***NEVER ACTUALLY DELETE*** this function.  It will serve as a good way to debug new sort methods
//vector<Point> Footprint3D::sort_unique_points(int cur_bin){
//	
//	//Find number of points in current bin and allocate space to store them in a working vector 'points'
//	int vec_size = (INTxx) current_points[cur_bin].size();
//	vector<Point> points;
//	points.assign(vec_size,Point());			//allocate the space needed as origin pts.  ALL will get overwritten.
//	
//	//If there are no points, return an empty vector
//	if (vec_size == 0) {
//		return points;	}
//	
//	//Load up the points that are in the current bin
//	for (int i = 0; i < vec_size; i++) {
//		points[i] = current_points[cur_bin][i];
////		if (points[i].get_y() < 1e-6) {
////			cout << "Shit... = " << points[i].get_y() << endl; }
//		points[i].set_id(i); 
////		cout << "point" << i << " = " << points[i] << endl;
//	
//	}
//	
//	//sort points by y value with a list
//	list<Point> mylist;
//	list<Point>::iterator it;
//	
//	bool need_origin = true;		//What's going on here?
//	int list_size;
//	it = mylist.begin();
//	
//	int first_index = 0;
//	
//	//Automatically accept the first point
//	mylist.insert(it, points[first_index]);			
//	
//	for (int i = first_index+1; i < vec_size; i++) {
//		
//		bool duplicate = false;
//		
//		list_size = (INTxx) mylist.size();
//		int j=0;	//makes sure we don't overflow the size of the present list
//		
//		double current_y = points[i].get_y();
//		double current_x = points[i].get_x();
//		
//		//check to see if the current pt is effectively a duplicate of another point already seen
//		it = mylist.begin();
//		for (int ix = 0; ix < list_size; ix++) {
//			double dist_between = sqrt(pow(current_y - it->get_y(),2) + pow(current_x - it->get_x(),2));
//			if (dist_between < 0.001) {
//				duplicate = true; 
//				break;
//			} 
//			it++; }
//		
//		if (!duplicate) {
//			//reset the iterator and find where the current point should fit in the list
//			it = mylist.begin();
//			while ((current_y <= it->get_y()) && (j < list_size)) { //finding where to place new point
//				
//				if (current_y == it->get_y()) {
//					//				cout << "TWO POINTS SHARE SAME YVAL!!!!!" << endl; 
//				}
//				
//				j++;
//				it++;
//			}
//		}
//		
//		
//		//BIG KLUDGE!!!! FIX THIS SOMETIME!!!
//		//having y=0 x>0 multiple points messes up algorithm, need to sort pts by xval for degenerate yvals
//		//until then, just zero out the x-parts for points with y=0
//		//overload comparison operators in point class.
//		
//		//This is still not totally right...sorting can still go awry near the end of the sorted list...
//		if ((points[i].get_y() == 0) && need_origin && !duplicate) {
//            cout << "I dont think this will ever realistically happen.  ERROR IF IT DOES THOUGH!!!!!!!!!~~~~~~~~~~~\n\n\n\n\n\n\n\n\n";
//            cout << points[i].get_y() << " " << need_origin << " " << duplicate << endl;
////			points[i].set_x(0); 
////			mylist.insert(it, points[i]);  //insert it into the list
////			need_origin = false;
//        }
//		else if ((list_size != 1) && (j+1 == list_size) && !duplicate) {	//put point at end of list
//			mylist.push_back(points[i]); }
//		else if (!duplicate) {				//Why do i care if it's positive????~~~~~~~~~~~~
//			mylist.insert(it, points[i]);  //insert it into the list (before element at 'it')
//		}
//		else {}
//	}
//	
//	it = mylist.begin();
//	vector<Point> remaining_points;
//	int num_remaining_pts = (INTxx) mylist.size();
//	remaining_points.assign(num_remaining_pts,Point());
//	
//	for (int i = 0; i < num_remaining_pts; i++) {
//		remaining_points[i] = *it;
//		it++; }
//	return remaining_points;
//	
//}











void Footprint3D::test_drive_swinging_arm(int timeIX) {
	// bin the points based on their z-location
	int num_current_points = (INTxx) all_points.size();

    // I think I need to allocate memory for the number of bins since this is known.  Why was this not an error earlier?
    current_points.assign(num_bins, vector<Point>());
	
	for (int i = 0; i < num_current_points; i++) {
		int bin_index = (INTxx) floor(all_points[i].get_z()/bin_size);
		if ((bin_index >=0) && (bin_index < num_bins)) {
//            cout << "current_points[" << bin_index << "].size = " << current_points[bin_index].size() << endl;
			current_points[bin_index].push_back(all_points[i]); }
		else {
			// If you wound up here, that means the altitude of the current point is out of the range you specified for the bins.
		} }
	
//	cout << "done binning by altitude" << endl;
	
	for (int cur_bin = 0; cur_bin < num_bins; cur_bin++) {
		vector<Point> unique_points = sort_unique_points(cur_bin);
        
//        cout << "up = [" << endl;
//        for (int up = 0; up < unique_points.size(); up++){
//            cout << unique_points[up].get_x() << "  " << unique_points[up].get_y() << ";" << endl;
//        }
//        cout << "];" << endl;
        
		vector<vector<int> > hull_indices = swing_the_arm(unique_points);	//pass by reference to avoid wasting time making a copy

		// Now need to store the footprint
		for (int cur_hull = 0; cur_hull < hull_indices.size(); cur_hull++) {
			footprint_storage3D[timeIX][cur_bin].push_back(vector<Point>());
			
			for (int ix = 0; ix < hull_indices[cur_hull].size(); ix++) {
				footprint_storage3D[timeIX][cur_bin][cur_hull].push_back(unique_points[hull_indices[cur_hull][ix]]); } }
	}
	
	
	return;
}



vector<vector<int> > Footprint3D::swing_the_arm(vector <Point> &unique_points) {

	vector<int> hull_indices;
	vector<vector<int> > hull_indices_storage;
	
	int vec_size = (INTxx) unique_points.size();
	
	vector<bool> available;
	available.assign(vec_size,true);	//make all pts available
	
	vector<int> avail_IX;				//available indices
	
	vector<bool> visited;
	visited.assign(vec_size, false);	//set all pts as unvisited
	
	int num_available;
	
	vector<Point> segments;
	vector<Point> temp_segments;
	
	// Why pushing back?  We know they're all available at this point...
	for (int i = 0; i < vec_size; i++) {
		if (available[i]) { 
			avail_IX.push_back(i); } } 
	num_available = (INTxx) avail_IX.size();
	
	int cur_hull = 0;
	while (num_available != 0) {
		hull_indices.clear();	//gets reset for every sub-hull
		hull_indices.push_back(avail_IX[0]);	//places highest y-point as first hull point (index)
		
		double angle_offset = 0;
		double last_angle = 0;
		
		int j=0;
		bool alive = true;
		
		
		int cur_pt;
		
		while (alive && ( j < 2*num_available ) ) {
			cur_pt = hull_indices.back();		//gets the most recent hull point (index)
			
			double last_best_angle = 0;
			double last_best_pt = cur_pt;
			
			//check to see which points are available
			//this kind of sucks, pull out the j==0 condition to before the loop, and do the else at the end of loop
			//	and make it a function
			
			
			double x1 = unique_points[cur_pt].get_x();
			double y1 = unique_points[cur_pt].get_y();
			
			//look through all available points for the one that's within arms length and 
			//  forms the smallest angle from the previous line segment
			for (int k = 0; k < num_available; k++) {
				double x2 = unique_points[avail_IX[k]].get_x();
				double y2 = unique_points[avail_IX[k]].get_y();
				
				double dx = x2-x1;
				double dy = y2-y1;
				
				double range = sqrt(dx*dx + dy*dy);
				
				if (cur_pt == avail_IX[k]) {	//throw out current point from calculations
					range = 1e6; }
				
				if ((range > 0.0) && (range < arm)) {
					//get angle and transform it be be from [0,2pi] measured clockwise from +x axis
					double angle = atan2(dy, dx);
					if (dy < 0) {
						angle = -angle; }
					else {
						angle = 2*PI - angle; }
					
					
					//some angles will appear to be smaller (better), but really they're nearly 360 degrees
					//	away from the last segment.  This catches those guys, as well as 360 deg rotations.
					//NOTE: MAY BE A BUG HERE!!! Only catches 360 rots on the upswing (offset neg) ???
					if (((fabs(angle)-eps) <= fabs(angle_offset)) && (angle_offset < 0)) {
						angle = angle + 2*PI; }
					
					double angle_from_previous_line = angle + angle_offset;
					if (((angle_from_previous_line) < last_best_angle) || (last_best_angle == 0)) {
						
						//add in the current segment
						temp_segments.push_back(unique_points[cur_pt]);
						temp_segments.push_back(unique_points[avail_IX[k]]);
						
						//do they intersect?
						int num_intersections = do_they_intersect(&temp_segments);
						
						//remove the current segment so we can check for others on the next iteration
						temp_segments.pop_back();
						temp_segments.pop_back();
						
						if ((num_intersections == 0) || ((num_intersections == 1) && (avail_IX[k] == hull_indices[0]))){
							last_best_angle = angle_from_previous_line;
							last_best_pt = avail_IX[k];
							last_angle = angle; }
						else { /*this isn't allowed, ignore it*/ } }
					
					
				} //ends range if
				
				
			} //ends k loop
			
			//save the good stuff
			segments.push_back(unique_points[hull_indices.back()]);
			segments.push_back(unique_points[last_best_pt]);
			hull_indices.push_back(last_best_pt);
			
			//this vector isnt actually necessary, but i think it makes me a little safer
			temp_segments.clear();
			temp_segments = segments;
			
			visited[hull_indices[j+1]] = true;
			
			angle_offset = PI - last_angle;
			
			if (hull_indices[0] == hull_indices.back()) {
				alive = false; }
			
			j++;
		}//ends while loop that creates a single footprint
		
		
		segments.clear();
		
		update_available_points(&available, &avail_IX, unique_points, hull_indices);
		num_available = (INTxx) avail_IX.size();
		
		//--------------------- Store the footprint as vector of hull_indices-------------
		
//		footprint_storage3D[cur_index][cur_bin].push_back(vector<Point>());
//		for (int i = 0; i < hull_indices.size(); i++) {
//			footprint_storage3D[cur_index][cur_bin][cur_hull].push_back(unique_points[hull_indices[i]]); }
		
		hull_indices_storage.push_back(vector<int>());
		for (int ix = 0; ix < hull_indices.size(); ix++) {
			hull_indices_storage[cur_hull].push_back(hull_indices[ix]); }

		cur_hull++;

	}//ends outer while loop that generates all footprints
	
//	cout << "swing: points in first hull = " << hull_indices_storage[0].size() << endl;
//	cout << "swing: num hulls = " << hull_indices_storage.size() << endl;

	return hull_indices_storage;
}



void Footprint3D::update_available_points(vector<bool> *available, vector<int> *avail_IX, vector<Point> &unique_points, vector<int> &hull_indices){
	//	cout << "in function, num_available = " << available->size() << endl;
	
	//******* This is only (atm) being designed to update after a single complete hull!!!!
	
	//	Point test = unique_points[hull_indices[0]];
	
	//find shape made by the completed footprint / closed region
	vector<Point> hull_segments;
	if (hull_indices.front() == hull_indices.back()) {
		int num_hull_pts = (INTxx) hull_indices.size();
		for (int i = 0; i < num_hull_pts-1; i++) {
			hull_segments.push_back(unique_points[hull_indices[i]]);
			hull_segments.push_back(unique_points[hull_indices[i+1]]);
		}	}
	else {
		cout << "Hull was not a closed shape.  Fatal Error." << endl;
		int num_hull_pts = (INTxx) hull_indices.size();
		cout << "badhull = [\n";
		for (int i = 0; i < num_hull_pts; i++) {
			cout << std::setprecision(12) << unique_points[hull_indices[i]].get_x() << "  "
                                            << unique_points[hull_indices[i]].get_y() << "  "
                                            << unique_points[hull_indices[i]].get_z() << endl; }
		cout << "];\n\n";
//		exit(982423);
	}
	
	//set the points that define the boundary as unavailable
	int hull_size = (INTxx) hull_indices.size();
	for (int i = 0; i < hull_size; i++) {
		available->at(hull_indices[i]) = false; }
	
	//i guess this is unnecessary since hull_segments is totally internal to this function, but whatever
	vector<Point> temp_segments = hull_segments;
	
	int num_available = (INTxx) avail_IX->size();	//hasnt been updated to include boundary pts, but doesn't matter
	for (int i = 0; i < num_available; i++) {
		Point reference_point(-1e4, unique_points[avail_IX->at(i)].get_y(),0,-5.);
		temp_segments.push_back(unique_points[avail_IX->at(i)]);
		temp_segments.push_back(reference_point);
		
		//check for intersections
		int num_intersections = do_they_intersect(&temp_segments);
		if ( num_intersections > 0) {
			if ((num_intersections % 2) == 1) {
				//point i is interior to the hull, don't want this one anymore
				available->at(avail_IX->at(i)) = false; }
			else {} }
		
		//remove the two points we just put on there
		temp_segments.pop_back();
		temp_segments.pop_back();
	}
	
	//this would be a good place to just update avail_IX since i have it anyways.
	vector<int> temp_avail = *avail_IX;		//copy the avail_IX vector (as a local object!)
	avail_IX->clear();						//then wipe it clean
	
	//cycle through points to find which are available and record the indices
	for (int i = 0; i < num_available; i++) {
		if (available->at(temp_avail[i])) {
			avail_IX->push_back(temp_avail[i]);	//populate it with the indices that are still available
		}
	}
	//clean up the temp vector
	temp_avail.clear();	
	
	
}


int Footprint3D::do_they_intersect(vector<Point> *segments) {
	/* Function designates the last segment as the important one and checks 
	 *	all previous segments against the final one; looking for an intersection.
	 *  All points are stored in *segments with the first two points defining seg1,
	 *  pts 3 and 4 make seg2, etc.
	 *	Return value is the number of times the final segment is intersected
	 */
	
	//how many segments?
	int num_pts = (INTxx) segments->size();
	
	int num_segs;
	if ((num_pts % 2) != 0) {
		cout << "ERROR!!!  You passed an incomplete segment!" << endl; 
		num_segs = -666;}
	else {
		num_segs = num_pts/2; }
	
	//	cout << "num_segs = " << num_segs << endl;
	//	cout << "first x val = " << segments->at(0).get_x() << endl;
	
	//going to overload the point class and use them like 2D vectors for a minute
	
	//checking against this segment (the last one)
	Point q = segments->at(2*(num_segs-1));
	Point s_prime = segments->at(2*(num_segs-1)+1);
	Point s = s_prime - q;
	
	//account for possible numerical errors
	double eps = 1e-12;
	
	bool printit = false;	//debugging

	int count = 0;
	for (int i = 0; i < num_segs-1; i++) {
		Point p = segments->at(2*i);
		Point r = segments->at(2*i+1) - p;
		

		if (printit) {
			cout << "segments{1} = [\n";
			for (int jx = 0; jx < num_pts-2; jx++) {
				cout << segments->at(jx).get_x() << "  " << segments->at(jx).get_y() << endl;
			}
			cout << "];\n\n";
			
			cout << "segments{2} = [\n";
			cout << segments->at(2*i).get_x() << "  " << segments->at(2*i).get_y() << endl;
			cout << segments->at(2*i+1).get_x() << "  " << segments->at(2*i+1).get_y() << endl;
			cout << "];\n\n";
			
			cout << "segments{3} = [\n";
			cout << segments->at(num_pts-2).get_x() << "  " << segments->at(num_pts-2).get_y() << endl;
			cout << segments->at(num_pts-1).get_x() << "  " << segments->at(num_pts-1).get_y() << endl;
			cout << "];\n\n";

			
		}
		
		double den =       r.get_x()*s.get_y() -     r.get_y()*s.get_x();
		double t_num = (q-p).get_x()*s.get_y() - (q-p).get_y()*s.get_x();
		double u_num = (q-p).get_x()*r.get_y() - (q-p).get_y()*r.get_x();
		
		double t = t_num/den;
		double u = u_num/den;
		
		//need to handle some special cases.  eps handles some numerical error
		//t = u = 0.  Share same beginning point
		if ((t < eps) && (u < eps)) {
			t = -50; }
		//t = 0, u = 1.  Points are end to end
		else if ((t < eps) && (fabs(u-1) < eps)) {
			t = -50; }
		//t = 1, u = 0.  Points are end to end
		else if ((u < eps) && (fabs(t-1) < eps)) {
			t = -50; }
		//t = 1, u = 1.  Share same end point
		else if ((fabs(t-1) < eps) && (fabs(u-1) < eps)) {
			t = -50; }
		//t = u = NaN.  Line segments are identical
		else if (isnan(t) && isnan(u)) {
			t = -666;
			u = -666; }
		//Something pathological happened, like one of the points exactly lies on the segment of interest
		else if ((t==1) || (t==0) || (u == 1) || (u == 0)){
			//Discard the current point
			count = 999;	//odd number will signal inside
			break; }
		//I think this means the segments overlap somehow.  This also shows up when one of the 'segments' has the same point for its beginning and end
		else if (isnan(1./fabs(t)) || isnan(1./fabs(u)) ) {
			//not sure what to do with this...treating it as not crossing
			t = -60;
			//cout << "segments overlap     overlapPoint = [" << segments->back().get_x() << " " << segments->back().get_y() << "];" << endl; 
		}
			
		//check to see if they intersect
		if ((t>=0) && (t<=1) && (u>=0) && (u<=1)) {
			count++; }	//they intersect!
		
	}
	
	
	return count ;
}

//void Footprint3D::append_to_existing_footprint(vector<vector<Point> > &total_points_at, double InitialUTC, double tstepMinutes,
//                                               double launchLat, double launchLon, double launchAzimuth){
//	Destruct_all_points();
//	
//	// Quick error check
//	if (allPointsFile.is_open()) {
//		cout << "ERROR: allPointsFile already open.  Please close the old points file manually to be sure this is what you want to do~~~~~~" << endl; 
//		//exit(232734);	
//	}
//	
//	//Making a copy even though not necessary ATM because there's no chance of the underlying vector falling out of scope...
//	//  This will, however, be safer in the future if I ever wind up doing a multithreaded code
//	all_points_passed_in = total_points_at;
//	
//	all_points_num_range = (INTxx) total_points_at.size();
//	all_points_UTC = InitialUTC;
//	all_points_delta_t = tstepMinutes;
//    all_points_launchLat = launchLat;
//    all_points_launchLon = launchLon;
//    all_points_launchAzimuth = launchAzimuth;
//	
//	append_to_existing_footprint();
//	return;
//}


void Footprint3D::append_to_existing_footprint(string newPointsFile){
	
	// load the new points to incorporate
	load_points_file(newPointsFile);		//Currently can only load one file at a time; no appending multiple files together
	
	append_to_existing_footprint();
	return;
}



//PRIVATE!!!
void Footprint3D::append_to_existing_footprint() {
	// We will use the existing footprint vector
	
	//NOTE: I had thought that maybe we could get away with only calling the swinging arm when there are new points to wrap up,
	//  but the issue is that you can't really tell when that's the case.  Checking intersections and finding that there are none
	//  only tells you that the boundary of the two shapes doesn't intersect...you could have two disconnected shapes.  So this
	//  function will have to run like a modified swinging arm where we sort the points, then instead of swinging right away do
	//  an update_available because we may find many interior points that we can kill to speed up the arm swing.
	

	// check that the delta_t's are the same	(MAYBE THIS SHOULD HAPPEN IN THE LOAD FUNCTION?!?!
	if (all_points_delta_t != footprint_delta_t) {
		cout << "ERROR: The points you're trying to append are incongruous with the existing footprint!~~~~~~~~~~~~" << endl;
		exit(233249); }
    
    // Check if launch info is the same
    if (footprint_launchLat == COMPOSITE_MISSION_LAT) {
        // Once a composite, always a composite.  Do nothing.
    } else if ((all_points_launchLat == footprint_launchLat) && (all_points_launchLon == footprint_launchLon)
               && (all_points_launchAzimuth == footprint_launchAzimuth)) {
        // They match, do nothing
    } else {
        // They don't match, so this is already a reentry footprint or is a possibly new composite
        if (footprint_launchLat == REENTRY_LAT){
            // Do nothing
        } else {
            // Composite mission
            footprint_launchLat = COMPOSITE_MISSION_LAT;
            footprint_launchLon = COMPOSITE_MISSION_LAT;
            footprint_launchAzimuth = COMPOSITE_MISSION_LAT;
            cout << "Youre starting a composite mission!" << endl;
        }
    }
    
    cout << "DEBUG YOU ARE HEREE!!!!!!!!!!!!" << endl;
	
	// Check the time range of the new points, may need to extend the time vector in either direction
	double all_points_final_UTC = all_points_UTC + all_points_delta_t*all_points_num_range/(3600*24);
	
	bool isLeading = (all_points_UTC < UTC_Initial);
	bool isTrailing = (all_points_final_UTC > UTC_Final);
	if (isLeading || isTrailing){
		//all_points is somehow outside the bounds of the existing footprint
		
		// how many leading time steps are needed as offset from existing first step
		int leadingTimeSteps = 0;
		if (isLeading) {
			leadingTimeSteps = (INTxx) ceil((UTC_Initial - all_points_UTC)*24*3600/footprint_delta_t);
			UTC_Initial = all_points_UTC; }
		
		// how many trailing time steps are needed as offset from existing final step
		int trailingTimeSteps = 0;
		if (isTrailing) { 
			//ceiling because no matter what, we need to add at least one timestep
			trailingTimeSteps = (INTxx) ceil((all_points_final_UTC - UTC_Final)*24*3600/footprint_delta_t);
			UTC_Final = all_points_final_UTC; }
		
		//---- Update the footprint vector
		
		//copy the existing footprint into a temp, delete the existing, reallocate with correct number of time steps
		vector<vector<vector<vector<Point> > > > tempFootprint;	        //[cur_index][z_bin][cur_hull][point]
		tempFootprint = footprint_storage3D;
		footprint_storage3D.clear();
		
		// Update the timing info (delta_t should be the same, leave it alone)
		footprint_num_range = leadingTimeSteps + ((INTxx) tempFootprint.size()) + trailingTimeSteps;
		footprint_UTC = UTC_Initial;
		
		cout << "lead = " << leadingTimeSteps << "  trail = " << trailingTimeSteps << "  fpnumrange = " << footprint_num_range << endl;
		
		// Allocate the updated footprint vector
		footprint_storage3D.assign(footprint_num_range, vector<vector<vector<Point> > >());
		for (int t = 0; t < footprint_num_range; t++) {
			footprint_storage3D[t].assign(num_bins,vector<vector<Point> >()); }
		
		// Copy the old info into the new vector at the appropriate time step
		for (int t = leadingTimeSteps; t < footprint_num_range - trailingTimeSteps; t++) {
			footprint_storage3D[t] = tempFootprint[t-leadingTimeSteps]; }
			
		// Delete the temp
		tempFootprint.clear();
        
	}
	
	
	
	
    cout << "Starting to append inside append function" << endl;
    // Set to 1000 ensures debug if-statements never hit
	static int debugCounter = 1000;
		
	// Loop over all the new all_points timesteps
	for (int all_points_tstep = 0; all_points_tstep < all_points_num_range; all_points_tstep++) {
		//
		destruct_leaves_of_current_points();
		
		
		load_points_at_timestep(all_points_tstep);
		
		// keep track of where there are new points
		vector <bool> newPointsHere(num_bins,false);
		
		
		// Bin the new points at the current tstep according to altitude
		for (int ix = 0; ix < all_points.size(); ix++) {
			int bin_index = (INTxx) floor(all_points[ix].get_z()/bin_size);
			if ((bin_index >=0) && (bin_index < num_bins)) {
				newPointsHere[bin_index] = true;
//                cout << "bin " << bin_index << "   ix " << ix << endl;
//                cout << "allpts[ix] = " << all_points[ix] << endl;
				current_points[bin_index].push_back(all_points[ix]); }
			else {
				// Throw this point out; it's outside the altitude range that we care about
			} }
						
		// The UTC that the current tstep corresponds to
		double currentUTC = all_points_UTC + all_points_delta_t*all_points_tstep/(3600*24); //xmarksthespot
		
		// The index in the updated footprint vector that the current UTC corresponds to
		int footprint_tstep = (INTxx) floor((currentUTC - UTC_Initial)*24*3600/footprint_delta_t);
				
		int debug_tstep = 40000*60;
		
		// Assuming points exist at this tstep / bin, find out how many points are exterior to the footprint
		for (int binIX = 0; binIX < num_bins; binIX++) {
			if (newPointsHere[binIX]) {
				
				//sorting only because I want to eliminate 'duplicate' points...don't actually need them in order yet
				vector<Point> unique_points = sort_unique_points(binIX);
				int num_hulls_here = (INTxx) footprint_storage3D[footprint_tstep][binIX].size();
				
				// DEBUG------- check to see what's really going on with matlab
				{
					if (footprint_tstep == debug_tstep) {
						cout << "%all_points_tstep = " << all_points_tstep << endl;
						
						for (int hullIX = 0; hullIX < num_hulls_here; hullIX++) {
							cout << "hull{" << debugCounter << "} = [\n";
							int num_points_here = (INTxx) footprint_storage3D[footprint_tstep][binIX][hullIX].size();
							
							for (int pt = 0; pt < num_points_here; pt++) {
								cout << std::setprecision(12);
								cout << "  " << footprint_storage3D[footprint_tstep][binIX][hullIX][pt].get_x()
								<< "  " << footprint_storage3D[footprint_tstep][binIX][hullIX][pt].get_y() << endl;
							}
							cout << "];\n\n";
						}
						
						if (debugCounter == 1) {
							cout << "points = [\n";					
							int current_points_here = (INTxx) unique_points.size();
							
							for (int ix = 0; ix < current_points_here; ix++) {
								cout << "  " << unique_points[ix].get_x() << "  " << unique_points[ix].get_y() << endl;
							}
							cout << "];\n\n";
						}
						
						debugCounter++;
					} //END DEBUG ------------------------
				}
				
				// The vectors that will be used to determine if new points are interior to existing footprint
				vector<vector<Point> > hullSegments;
				hullSegments.assign(num_hulls_here,vector<Point>());
				
				// The reference point
				Point referencePt(-25000., -25000, 0, -5.);
				

				for (int hullIX = 0; hullIX < num_hulls_here; hullIX++) {
					int num_points_here = (INTxx) footprint_storage3D[footprint_tstep][binIX][hullIX].size();	//NOTE: closed shape --> last point = first point
					
					// Note that the last two points of hullSegments are for the point we're trying to determine if enclosed
					for (int pt = 0; pt < num_points_here-1; pt++) {
                        hullSegments[hullIX].push_back(footprint_storage3D[footprint_tstep][binIX][hullIX][pt]);
						hullSegments[hullIX].push_back(footprint_storage3D[footprint_tstep][binIX][hullIX][pt+1]);
					}
					
					// Check to make sure it's a closed shape
					if (hullSegments[hullIX][0] == hullSegments[hullIX][2*(num_points_here-2)+1]) {
						//cout << "It's a closed shape!!!" << endl; 
					}
					else {
						cout << "BOOOOOOO: Footprint not a closed shape!!!" << endl; }
					

					// Push back the reference point
					hullSegments[hullIX].push_back(referencePt);

				}//Ends hullIX loop

				// Push the points to check into hullSegments and see if there are intersections
				int unique_points_here = (INTxx) unique_points.size();
				vector<Point> externalPts;
				
				for (int newpt = 0; newpt < unique_points_here; newpt++) {
					
					// Assume point is external unless we find evidence to the contrary
					bool isExternal = true;	
					
					for (int hullIX = 0; hullIX < num_hulls_here; hullIX++) {
						// Push back the current point so we can test it
						hullSegments[hullIX].push_back(unique_points[newpt]);

						//check for intersections
						int num_intersections = do_they_intersect(&hullSegments[hullIX]);
						
						if ((num_intersections % 2) == 1) {
							//point i is interior to the hull, throw it out
							isExternal = false;
							
							//remove the point we just put on there so we can break safely
							hullSegments[hullIX].pop_back();
							
							//don't need to check any more hulls, we know this point is internal
							break;
						}
						else {
							//this point is exterior to the hull, keep it (i.e. do nothing)
							
							//remove the point we just put on there, doing this here so the break above works cleanly
							Point pointOfInterest = hullSegments[hullIX].back();
							hullSegments[hullIX].pop_back(); 
							
							// Want to catch points that are effectively hull points
							int num_points_existing = (INTxx) footprint_storage3D[footprint_tstep][binIX][hullIX].size();	//NOTE: closed shape --> last point = first point
						
							for (int existingPt = 0; existingPt < num_points_existing; existingPt++) {
								Point currExisting = footprint_storage3D[footprint_tstep][binIX][hullIX][existingPt];
								double dist = sqrt(pow(pointOfInterest.get_x() - currExisting.get_x(),2) + pow(pointOfInterest.get_y() - currExisting.get_y(),2));
								if (dist < 0.0001) {
									//These points are the same! count the point as internal so we can throw it out
									isExternal = false; } }
							
							// Want to catch when we have existing footprints of only a single point
						
						
						} 
					}
					
					// Keep the point if we've determined it to be external
					if (isExternal) {
						externalPts.push_back(unique_points[newpt]); }
						
				}//Ends for newpt loop
					
                // NOTE: It may be worth breaking the function up here and writing the new external points to a file or storing in a vec
                //  if we're confident that there will be very few new points to consider.  This way we could wrap all the new points
                //  up at once, later, rather than swinging an arm every time we find a SINGLE new point.
					
				// NOTE: It occurs to me that it shouldn't matter if the points are sorted into order or not...you only need to make
				//  sure that the point you start with is itself on the boundary.  So you could just look for the first point and then
				//  let the others be where they fell.  If there is another hull, then simply look through the remaining points for 
				//  the next important point...
				
				// NOTE:  It appears that if you run the same points twice, it will think that most of the boundary pts are external.
				//  I'm okay with this for now, but soon that should be fixed.  Can just check the external pts for distance from hull pts.
					
				
				if ((footprint_tstep == debug_tstep) && (debugCounter <= 3)){
					//First time through counter will be 2 because it already got incremented earlier
					cout << "Extpoints = [\n";					
					int external_points_here = (INTxx) externalPts.size();
					
					for (int ix = 0; ix < external_points_here; ix++) {
						cout << "  " << externalPts[ix].get_x() << "  " << externalPts[ix].get_y() << endl;
					}
					cout << "];\n\n";
				}
                
				if (externalPts.size() > 0) {					
					//load all of the points into current_points so we can sort them again
					current_points[binIX].clear();
					current_points[binIX] = externalPts;
					
					//Collect the footprint points and the new interior points and swing the arm
					for (int hullIX = 0; hullIX < num_hulls_here; hullIX++) {
						for (int existingPt = 0; existingPt < footprint_storage3D[footprint_tstep][binIX][hullIX].size(); existingPt++) {
							//load all of the points into current_points so we can sort them again
							current_points[binIX].push_back(footprint_storage3D[footprint_tstep][binIX][hullIX][existingPt]); } }
														
					if ((footprint_tstep == debug_tstep) && (debugCounter <= 3)){
						//First time through counter will be 2 because it already got incremented earlier
						cout << "CurrentPoints = [\n";											
						for (int ix = 0; ix < current_points[binIX].size(); ix++) {
							cout << "  " << current_points[binIX][ix].get_x() << "  " << current_points[binIX][ix].get_y() << endl;
						}
						cout << "];\n\n";
					}
					
					//sort them again and update the unique points vector
					unique_points = sort_unique_points(binIX);
					
					if ((footprint_tstep == debug_tstep) && (debugCounter <= 3)){
						//First time through counter will be 2 because it already got incremented earlier
						cout << "UniquePoints = [\n";											
						for (int ix = 0; ix < unique_points.size(); ix++) {
							cout << "  " << unique_points[ix].get_x() << "  " << unique_points[ix].get_y() << endl;
						}
						cout << "];\n\n";
					}
					
					//Swing it
					vector<vector<int> > hull_indices = swing_the_arm(unique_points);	//pass by reference to avoid wasting time making a copy
					//cout << "arm's been swung" << endl;
					
					// Clear out the old footprint at this point
					footprint_storage3D[footprint_tstep][binIX].clear();
					
					// Now need to store the updated footprint
					for (int cur_hull = 0; cur_hull < hull_indices.size(); cur_hull++) {
						footprint_storage3D[footprint_tstep][binIX].push_back(vector<Point>());
						
						for (int ix = 0; ix < hull_indices[cur_hull].size(); ix++) {
							footprint_storage3D[footprint_tstep][binIX][cur_hull].push_back(unique_points[hull_indices[cur_hull][ix]]); } }
					
				}//ends if(ext.size() > 0)
					
				
					
					
					// Finish and clean up
					hullSegments.clear();
					externalPts.clear();
				
			} //ends if(newpointshere)
			else {
				//cout << "No points here!" << endl; 
			}
		}//ends bin loop
	}//ends all_points_tstep loop
	
    cout << "DONE with append inside append function" << endl;

	return;
}














/*!
 * This gets use for two things.  First, you can take an existing footprint
 *      and change the timestep size to be bigger by merging FPs that fall within
 *      the new deltaT.  A new footprint is generated from the merged ones.
 *      Second, if you've done some merging yourself, outside of this function,
 *      and want to remove possibly duplicate points or find the new outter 
 *      boundaries of the shape, you can call this without any deltaT change
 *      and it will rerun the footprint generation algorithm for you.
 */
void Footprint3D::SmoothedOut(double newDeltaT){

    // ------------
    // This first big loop simply takes the incoming footprint structure, and converts it into Points
    // ------------
    
    // This is where we'll temporarily put the points
    vector<vector<Point> > all_points_temp;	//[timestep][point index]
    
    // Loop through all the points and bin them only by time INDEX (int)
    int numTimeStepsFP = (INTxx) footprint_storage3D.size();
    all_points_temp.assign(numTimeStepsFP, vector<Point>());
    for (int tx = 0; tx < numTimeStepsFP; tx++){
        int curPtsHere = 0;
        
        // Find out how many points there are at this timestep (basically just counting so we can allocate the vector)
        int numBins = (INTxx) footprint_storage3D[tx].size();
        for (int Bx = 0; Bx < numBins; Bx++){
            
            int numHulls = (INTxx) footprint_storage3D[tx][Bx].size();
            for (int Hullx = 0; Hullx < numHulls; Hullx++){
                curPtsHere += (INTxx) footprint_storage3D[tx][Bx][Hullx].size(); } }
        
        // Allocate the space
        Point EmptyPoint;
        all_points_temp[tx].assign(curPtsHere, EmptyPoint);
        
        // Load up the points
        int curIndex = 0;
        for (int Bx = 0; Bx < numBins; Bx++){
            int numHulls = (INTxx) footprint_storage3D[tx][Bx].size();
            
            for (int Hullx = 0; Hullx < numHulls; Hullx++){
                int numPts = (INTxx) footprint_storage3D[tx][Bx][Hullx].size();
                
                for (int Px = 0; Px < numPts; Px++){
                    all_points_temp[tx][curIndex] = footprint_storage3D[tx][Bx][Hullx][Px];
                    curIndex++;     } } }
    }   // ends tx loop
    
//    cout << "c++ : " << numTimeStepsFP << endl;
    
    // ------------
    // Now, if the timestep has changed, we rebin the points according to that timestep.  Otherwise, we do nothing.
    // ------------
    
    // If the incoming deltaT is greater, then not only will we smooth the points, but we'll also adjust the time binning
    vector<vector<Point> > all_points_REBIN;	//[timestep][point index]
    int numTimeSteps = 0;
    
    if (newDeltaT > footprint_delta_t){
        int numTimeStepsAllPoints = ceil(numTimeStepsFP * footprint_delta_t/newDeltaT);
        all_points_REBIN.assign(numTimeStepsAllPoints, vector<Point>());
        
//        cout << "c++ : numTimeStepsFP " << numTimeStepsFP << endl;
//        cout << "c++ : footprint_delta_t " << footprint_delta_t << endl;
//        cout << "c++ : newDeltaT " << newDeltaT << endl;
//        cout << "c++ : numTimeStepsAllPoints " << numTimeStepsAllPoints << endl;

        // Loop through all the points that you just loaded up
        for (int tx = 0; tx < numTimeStepsFP; tx++){
            int newTX = floor(tx*footprint_delta_t/newDeltaT);
//            cout << "c++ :      newTX " << newTX << endl;

            
//            // Push back the vector of points from this original tx
//            all_points_REBIN[newTX].push_back(all_points_temp[tx]);
            
            // Push back is for single elements, thus not valid here.  Need to insert
            all_points_REBIN[newTX].insert(all_points_REBIN[newTX].end(), all_points_temp[tx].begin(), all_points_temp[tx].end());
        }
        
        // Pass these off
        numTimeSteps = numTimeStepsAllPoints;
        footprint_num_range = numTimeSteps;
        footprint_delta_t = newDeltaT;
        all_points_passed_in = all_points_REBIN; }
    else {
        
        // No timestep change, use the points you already found
        numTimeSteps = numTimeStepsFP;
        all_points_passed_in = all_points_temp;
    }
    
//    cout << "c++ : " << numTimeSteps << endl;

    // ------------
    // Finally, take the points and rerun the footprint algorithm.
    // ------------
    
    // Clear the old footprint
    footprint_storage3D.clear();
    
    // Allocate the footprint vector
	footprint_storage3D.assign(numTimeSteps,vector<vector<vector<Point> > >());
	for (int t = 0; t < footprint_num_range; t++) {
		footprint_storage3D[t].assign(num_bins,vector<vector<Point> >()); }
    
    // Now fill it up
    generate_footprint_at_timesteps();
    
    return;
}






// Cython deals in pointers, not the actual objects themselves
void Footprint3D::MergeFootprintVectors(Footprint3D *incomingFP){
//    cout << "hopefully this works" << endl;
    MergeFootprintVectors(*incomingFP);
}


void Footprint3D::MergeFootprintVectors(Footprint3D &incomingFP){
    
            //	// -------- Copy the vitals of the existing footprint (call it FP2)
            //	int FP2_footprint_num_range = footprint_num_range;
            //	double FP2_footprint_UTC = footprint_UTC;
            //	double FP2_footprint_delta_t = footprint_delta_t;
            //	double FP2_bin_size = bin_size;
            //	int FP2_num_bins = num_bins;
            //    
            //    //	// Not vital, but will want to use it.
            //    //	double FP2_UTC_Final = UTC_Final;
            //	double FP2_finalUTC = FP2_footprint_UTC + FP2_footprint_num_range*FP2_footprint_delta_t/(24.*60.);
            //    
            //	vector<vector<vector<vector<Point> > > > FP2_storage3D;	//[cur_index][z_bin][cur_hull][point]
            //	FP2_storage3D = footprint_storage3D;

    
            //	// --------- Wipe out the old footprint and load up the new one
            //	load_footprint_as_vector(new_footprint_name);
            //    
            //	// Calculate a few things that the main constructor calculates just to be consistent
            //	double UTC_Initial = footprint_UTC;
            //	double delta_t = footprint_delta_t;
            //	double UTC_Final = UTC_Initial + delta_t*footprint_num_range/(60*24);
            //    
            //	// What's the total range?
            //	int TEMP_num_range = (INTxx) round(fabs(UTC_Final - FP2_finalUTC)*24*60/delta_t);
            //	cout << "tempnumrange = " << TEMP_num_range << endl;
    
    // ******
    // The above should all not be necessary because we're passing in the footprint to merge as its own object
    // ******
    
    // MUST BE CAREFUL with translation.  Anything that refers to FP2 gets mapped to the self (cur) footprint.
    //                                    Anything that refers to the old self footprint gets mapped to incomingFP
    
    double cur_UtcFinal = this->getFinalUTC();// UTC_Final;
    double in_UtcFinal = incomingFP.getFinalUTC();
    
    double cur_NumRange = this->getNumRange();
    double in_NumRange = incomingFP.getNumRange();
    
    double cur_InitialUtc = this->getInitialUTC();
    double in_InitialUtc = incomingFP.getInitialUTC();
    
    double cur_DeltaT = this->getDeltaT();      //seconds
    if (cur_DeltaT != incomingFP.getDeltaT()) {
        cout << "The delta_t's are not equal, thus cannot merge.  Exiting function without merging.  ERROR!!!!" << endl;
        return; }
    
    double cur_NumBins = this->getNumBins();
    if (cur_NumBins != incomingFP.getNumBins()) {
        cout << "The number of bins are not equal, thus cannot merge.  Exiting function without merging.  ERROR!!!!" << endl;
        return; }
    
    double cur_BinSize = this->getBinSize();
    if (cur_BinSize != incomingFP.getBinSize()) {
        cout << "The size of the bins are not equal, thus cannot merge.  Exiting function without merging.  ERROR!!!!" << endl;
        return; }
    
    double TEMP_num_bins = cur_NumBins;     // Leave the Temp_ because later we may want to allow for non-equal bins to be merged
    
    double cur_launchLat = this->getLaunchLat();
    double in_launchLat = incomingFP.getLaunchLat();
    
    double cur_launchLon = this->getLaunchLon();
    double in_launchLon = incomingFP.getLaunchLon();
    
    double cur_launchAzimuth = this->getLaunchAzimuth();
    double in_launchAzimuth = incomingFP.getLaunchAzimuth();
    
//    cout << "cur_launchLat = " << this->getLaunchLat() << endl;
//    cout << "in_launchLat = " << in_launchLat << endl;
//    cout << "cur_launchLon = " << this->getLaunchLon() << endl;
//    cout << "in_launchLon = " << in_launchLon << endl;
//    cout << "cur_launchAzimuth = " << this->getLaunchAzimuth() << endl;
//    cout << "in_launchAzimuth = " << in_launchAzimuth << endl;
    
    double temp_launchLat = cur_launchLat;
    double temp_launchLon = cur_launchLon;
    double temp_launchAzimuth = cur_launchAzimuth;
    
    // Check if launch info is the same
    if (cur_launchLat == COMPOSITE_MISSION_LAT) {
        // Do nothing, once you're composite there's no going back.
//    } else if ((cur_launchLat == in_launchLat) && (cur_launchLon == in_launchLon) && (cur_launchAzimuth == in_launchAzimuth)) {
    } else if ((cur_launchLat == in_launchLat) && (cur_launchLon == in_launchLon)) {
        // They match, do nothing
    } else { 
        // They don't match, so this is already a reentry footprint or is a possibly new composite
        if (cur_launchLat == REENTRY_LAT){
            // Do nothing for now
        } else {
            
            cout << "DEBUG" << endl;
            cout << cur_launchLat << "   " << in_launchLat << endl;
            cout << cur_launchLon << "   " << in_launchLon << endl;
            cout << cur_launchAzimuth << "   " << in_launchAzimuth << endl;
            
            // Composite mission
            temp_launchLat = COMPOSITE_MISSION_LAT;
            temp_launchLon = COMPOSITE_MISSION_LAT;
            temp_launchAzimuth = COMPOSITE_MISSION_LAT;
            cout << "Youre starting a composite mission!" << endl;
        }
    }
    
    

	// What's the total range?
//	if (in_UtcFinal > cur_UtcFinal) {
//		TEMP_num_range = in_NumRange; }
//	else {
//		TEMP_num_range = cur_NumRange; }
    
    // These values will get updated if necessary in the (isLeading || isTrailing) logic
    double TEMP_InitialUtc = cur_InitialUtc;
    double TEMP_FinalUtc = cur_UtcFinal;
    
    int leadingTimeSteps = 0;
    int trailingTimeSteps = 0;
    
    int FP1start = 0;
	int FP2start = 0;
    
    bool isLeading = (in_InitialUtc < cur_InitialUtc);
	bool isTrailing = (in_UtcFinal > cur_UtcFinal);
    
    if (isLeading) {
        leadingTimeSteps = (INTxx) round((cur_InitialUtc - in_InitialUtc)*24*3600/cur_DeltaT);
        TEMP_InitialUtc = in_InitialUtc;
        FP1start = 0;
        FP2start = leadingTimeSteps; }
    else {
        FP2start = 0;
        FP1start = (INTxx) round((-cur_InitialUtc + in_InitialUtc)*24*3600/cur_DeltaT); }
    
    
    if (isTrailing) {
			//ceiling because no matter what, we need to add at least one timestep
        trailingTimeSteps = (INTxx) ceil((in_UtcFinal - cur_UtcFinal)*24*3600/cur_DeltaT);
        TEMP_FinalUtc = in_UtcFinal; }
    
    
//    cout << "in merge, leading = " << leadingTimeSteps << "   trailing = " << trailingTimeSteps << endl;
    
//	if (isLeading || isTrailing){
//		//all_points is somehow outside the bounds of the existing footprint
//		
//		// how many leading time steps are needed as offset from existing first step
//		if (isLeading) {
////			leadingTimeSteps = (INTxx) ceil((cur_InitialUtc - in_InitialUtc)*24*60/cur_DeltaT);
//			leadingTimeSteps = (INTxx) round((cur_InitialUtc - in_InitialUtc)*24*60/cur_DeltaT);
//			TEMP_InitialUtc = in_InitialUtc; }
//		
//		// how many trailing time steps are needed as offset from existing final step
//		if (isTrailing) {
//			//ceiling because no matter what, we need to add at least one timestep
//			trailingTimeSteps = (INTxx) ceil((in_UtcFinal - cur_UtcFinal)*24*60/cur_DeltaT);
//			TEMP_FinalUtc = in_UtcFinal; } }
    
    int TEMP_num_range = leadingTimeSteps + cur_NumRange + trailingTimeSteps;
//    footprint_num_range = leadingTimeSteps + ((INTxx) tempFootprint.size()) + trailingTimeSteps;

    
	// ---------- Now start to combine them into a single temp_vector



//	if (in_InitialUtc == cur_InitialUtc){
//		cout << "The UTC's are equal!" << endl; }
//	else if (in_InitialUtc > cur_InitialUtc) {
//		// This means that FP2 starts first
//		FP2start = 0;
//        
//		// How many delta_t's before FP1 starts?
//		cout << "fp1start = " << (in_InitialUtc - cur_InitialUtc)*24*60 / cur_DeltaT << endl;
//		FP1start = (INTxx) round((in_InitialUtc - cur_InitialUtc)*24*60 / cur_DeltaT);
//		cout << "again but as an int = " << FP1start << endl;
//        
//	}
//	else if (in_InitialUtc < cur_InitialUtc) {
//		// This means that FP1 starts first
//		FP1start = 0;
//        f
//		// How many delta_t's before FP2 starts?
//		cout << "fp2start = " << (cur_InitialUtc - in_InitialUtc)*24*60 / cur_DeltaT << endl;
//		FP2start = (INTxx) round((cur_InitialUtc - in_InitialUtc)*24*60 / cur_DeltaT);
//		cout << "again but as an int = " << FP2start << endl;
//	}
//	else {
//		cout << "The UTC's for the two footprints aren't equal, this case isn't handled yet";
//		exit(9); }


	vector<vector<vector<vector<Point> > > > TEMP_storage3D;	//[cur_index][z_bin][cur_hull][point]
	TEMP_storage3D.assign(TEMP_num_range, vector<vector<vector<Point> > >());
    
    
	// Do the merge
	for (int t = 0; t < TEMP_num_range; t++) {
		TEMP_storage3D[t].assign(TEMP_num_bins, vector<vector<Point> >());
        
		// These indices will keep track of the separate footprints
		int t1 = t - FP1start;  //incomingFP
		int t2 = t - FP2start;  //this* 

		for (int binx = 0; binx < cur_NumBins; binx++) {
			int existing_hulls = 0;
			int merging_hulls = 0;
            
			if ((t1 >= 0) && (t1 < in_NumRange)) {
                merging_hulls= (INTxx) incomingFP.footprint_storage3D[t1][binx].size(); }
            
			if ((t2 >= 0) && (t2 < cur_NumRange)){
                existing_hulls = (INTxx) this->footprint_storage3D[t2][binx].size(); }

			int num_hulls = merging_hulls + existing_hulls;
			TEMP_storage3D[t][binx].assign(num_hulls, vector<Point>());
            
			int curr_hull = 0;
			for (int hullx = 0; hullx < merging_hulls; hullx++) {
				int num_points = (INTxx) incomingFP.footprint_storage3D[t1][binx][hullx].size();
//				cout << "num poins = " << num_points << endl;
				TEMP_storage3D[t][binx][curr_hull].assign(num_points, Point());
				memcpy((void *) &TEMP_storage3D[t][binx][curr_hull][0], (void *) &incomingFP.footprint_storage3D[t1][binx][hullx][0], num_points*sizeof(Point));
				curr_hull++; }
            
			for (int hullx = 0; hullx < existing_hulls; hullx++) {
				int num_points = (INTxx) this->footprint_storage3D[t2][binx][hullx].size();
//				cout << "num pointtts = " << num_points << endl;
				TEMP_storage3D[t][binx][curr_hull].assign(num_points, Point());
				memcpy((void *) &TEMP_storage3D[t][binx][curr_hull][0], (void *) &this->footprint_storage3D[t2][binx][hullx][0], num_points*sizeof(Point));
				curr_hull++; }
		}
	}
    
	// ----- Wipe out existing footprint and replace with new one
	Destruct_footprint();
    
    // <-- Left off here...should really define some local variables instead!
//    vector<vector<Point> > assemble_all_points_debris(double tstepMinutes);


	// -------- Copy the vitals of the footprint
	footprint_num_range = TEMP_num_range;
	footprint_UTC = TEMP_InitialUtc;
	footprint_delta_t = cur_DeltaT; //seconds
	bin_size = cur_BinSize;
	num_bins = cur_NumBins;
    
    footprint_launchLat = temp_launchLat;
    footprint_launchLon = temp_launchLon;
    footprint_launchAzimuth = temp_launchAzimuth;
    
	footprint_storage3D = TEMP_storage3D;
    
	// Calculate a few things that the main constructor calculates just to be consistent
	UTC_Initial = footprint_UTC;
//	delta_t = cur_DeltaT;
	UTC_Final = UTC_Initial + footprint_delta_t*TEMP_num_range/(3600*24);
    
//    cout << "cur_launchLat000 = " << this->getLaunchLat() << endl;

    
	// Write the merged vector to a file (if asked to)
//	store_footprint_as_vector(FPOut_name);
    
//	cout << "leaving merge" << endl;
    
    
    return;
}

int Footprint3D::ProjectAllPointsDown(){
    
//    cout << "launcLatTEMP = " << getLaunchLat() << endl;

    
    // Be safe and clear this out
    all_points.clear();
    all_points = CollapseFootprintToPoints();      // Compress all the points down to one timestep
    
    // Clear out the existing footprint
    footprint_storage3D.clear();    //i think since this is a std::vector, this should be safe enough
    
    int old_num_range = footprint_num_range;
    
    // Set the relevant all_points parameters
    all_points_num_range = 1;
//    all_points_delta_t = footprint_delta_t*footprint_num_range;
    all_points_delta_t = footprint_delta_t;  //seconds
    
    bin_size = 20.; //set it above the NAS to capture everything
    num_bins = 1;
    
    // Set those to be true for the footprint
    footprint_num_range = all_points_num_range;
    footprint_delta_t = all_points_delta_t;     //seconds
    
    // re-Allocate the footprint vector
	footprint_storage3D.assign(footprint_num_range,vector<vector<vector<Point> > >());
	for (int t = 0; t < footprint_num_range; t++) {
		footprint_storage3D[t].assign(num_bins,vector<vector<Point> >()); }
    
    test_drive_swinging_arm(0); // zero because we just compressed everything down to one tstep
    
    
    
    // Set the delta_t to encompass all the points
    // Since we haven't generated a footprint yet, the points file will dictate the timing info of the footprint we're about to create
//	footprint_num_range = all_points_num_range;
//	footprint_UTC = all_points_UTC;
//	footprint_delta_t = all_points_delta_t;
    
    return old_num_range;
}

// Preserves the time information
vector<vector<Point> > Footprint3D::DumpFootprintToPoints(){
	// Get the points into the appropriate vector
	std::vector<std::vector<Point> > total_points_at;	//trajectory and debris points all put together by timestep
	total_points_at.assign(footprint_num_range,std::vector<Point>());

	for (int t = 0; t < footprint_num_range; t++){
		for (int ix = 0; ix < footprint_storage3D[t].size(); ix++) {
			for (int jx = 0; jx < footprint_storage3D[t][ix].size(); jx++) {
				for (int kx = 0; kx < footprint_storage3D[t][ix][jx].size(); kx++) {
					total_points_at[t].push_back(footprint_storage3D[t][ix][jx][kx]); }}}}
    
    return total_points_at;
}

// Dumps all the points to a single timestep.  
vector<Point> Footprint3D::CollapseFootprintToPoints(){
	// Get the points into the appropriate vector
	std::vector<Point> total_points_at;	//trajectory and debris points all put together by timestep
//	total_points_at.assign(1,std::vector<Point>());
    
	for (int t = 0; t < footprint_num_range; t++){
		for (int ix = 0; ix < footprint_storage3D[t].size(); ix++) {
			for (int jx = 0; jx < footprint_storage3D[t][ix].size(); jx++) {
				for (int kx = 0; kx < footprint_storage3D[t][ix][jx].size(); kx++) {
					total_points_at.push_back(footprint_storage3D[t][ix][jx][kx]); }}}}
    
    return total_points_at;
}




void Footprint3D::AddToFootprintUTC(double plus_hours, double plus_minutes, string outputFileName) {
	// This is NOT the way this should be handled, but I don't have the time to fix the issue and debug it
	// SHOULD totally drop the UTC thing and store information in local time.

	UTC_Initial = UTC_Initial + plus_hours/24. + plus_minutes/(60*24);
	UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);

	// Have the option of writing the new footprint to a file
	if (outputFileName != "none") {
		store_footprint_as_vector(outputFileName); }

	return;
}

void Footprint3D::ChangeLaunchSiteToDeg(double gdlatIN, double lonIN){
    ChangeLaunchSiteToRad(gdlatIN * (PI/180), lonIN * (PI/180));
    return;
}


void Footprint3D::ChangeLaunchSiteToRad(double gdlatIN, double lonIN){
    
    Point oldPoint(footprint_launchLat * (PI/180), footprint_launchLon * (PI/180), 0, 0);
    Point deltaPoint(gdlatIN, lonIN, 0, 0);
    
//    cout << "oldPoint = " << footprint_launchLat << ", " << footprint_launchLon << endl;
    deltaPoint = deltaPoint - oldPoint;

    //	vector<vector<vector<vector<Point> > > > footprint_storage3D;	//[cur_index][z_bin][cur_hull][point]
    for (int curTime = 0; curTime < footprint_num_range; curTime++){
        for (int curZ = 0; curZ < num_bins; curZ++){
            
            int numHullsHere = footprint_storage3D[curTime][curZ].size();
            for (int curHull = 0; curHull < numHullsHere; curHull++){
                
                int numPtsHere = footprint_storage3D[curTime][curZ][curHull].size();
                for (int curPt = 0; curPt < numPtsHere; curPt++){

                    footprint_storage3D[curTime][curZ][curHull][curPt]
                    = footprint_storage3D[curTime][curZ][curHull][curPt] + deltaPoint;  } } } }
    
    // Make sure to update the launch coords!
    Point newLaunchLocation = oldPoint + deltaPoint;
    footprint_launchLat = newLaunchLocation.get_gdLatDeg();
    footprint_launchLon = newLaunchLocation.get_lonDeg();
    
//    cout << "new launch lat and lon = " << footprint_launchLat << "   " << footprint_launchLon << endl;
//    cout << "was shooting for       = " << gdlatIN*(180/PI)  << "    " << lonIN*(180/PI) << endl;
    
    return;
}
 
void Footprint3D::ChangeAzimuthToDeg(double newAzimuth){
    ChangeAzimuthToRad(newAzimuth * (PI/180));
    return;
}

void Footprint3D::ChangeAzimuthToRad(double newAzimuth){
    
    Point refPoint(footprint_launchLat * (PI/180), footprint_launchLon * (PI/180), 0, 0);
    double rotAngle = -(newAzimuth - footprint_launchAzimuth*(PI/180));  //Negative sign because Azimuth is measured clockwise from North
                                                                //  but normal angle convention is to measure counterclockwise.
    cout << "rotation Angle = " << rotAngle * (180/PI) << endl;

    //	vector<vector<vector<vector<Point> > > > footprint_storage3D;	//[cur_index][z_bin][cur_hull][point]
    for (int curTime = 0; curTime < footprint_num_range; curTime++){
        for (int curZ = 0; curZ < num_bins; curZ++){
            
            int numHullsHere = footprint_storage3D[curTime][curZ].size();
            for (int curHull = 0; curHull < numHullsHere; curHull++){
                
                int numPtsHere = footprint_storage3D[curTime][curZ][curHull].size();
                for (int curPt = 0; curPt < numPtsHere; curPt++){
                    
                    footprint_storage3D[curTime][curZ][curHull][curPt].rotateThisAboutAnotherPoint(refPoint, rotAngle);
                      } } } }
    
    // Make sure to update the new launch azimuth!
    footprint_launchAzimuth = newAzimuth*(180./PI); 
    
    return;
}


double Footprint3D::ChopTimeAt(double seconds){
    
    int numChoppedSteps = ceil(seconds / footprint_delta_t);    // That is, number of steps left AFTER chopping
    cout << "numChoppedSteps = " << numChoppedSteps << endl;
    cout << "footprint_num_range = " << footprint_num_range << endl;
    cout << "footprint_delta_t (in seconds) = " << footprint_delta_t << endl;
    
    if (numChoppedSteps < footprint_num_range){
        footprint_num_range = numChoppedSteps;
        footprint_storage3D.resize(footprint_num_range);
        
        //UTC_Intial stays the same, but final changes
        UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);
                
    } else if (footprint_num_range == 1) {
        cout << "WARNING!!!  Only see a single timestep.  *** Resizing that delta_t ***" << endl;
        footprint_delta_t = ceil(seconds);
        
        //UTC_Intial stays the same, but final changes
        UTC_Final = UTC_Initial + footprint_delta_t*footprint_num_range/(3600*24);
        
    } else {
        cout << "WARNING!!!  You're trying to chop some seconds but it's too many.  Doing Nothing." << endl;
    }
    
    // Return the time length of the current footprint
    return footprint_num_range*footprint_delta_t;
}



void Footprint3D::SlideFootprintBySeconds(int howManySeconds){
    // Only going to allow forward slides for now
    if (howManySeconds < 0){
        cout << "\n\nERROR: ONLY ALLOWED TO SLIDE FORWARD. Exiting\n";
        exit(5);    }
    
    // Make a copy of footprintStorage
//    Footprint3D oldFootprint_storage3D;
//    oldFootprint_storage3D              = *this;
    vector<vector<vector<vector<Point> > > > oldFootprint_storage3D = footprint_storage3D;
    int oldNumRange                     = footprint_num_range;
    double oldDeltaT                    = footprint_delta_t;
    
    // Make a new footprintStorage that has the number of timesteps that we want
    int addedSteps                      = ceil(howManySeconds/oldDeltaT);
    int desiredNumRange                 = oldNumRange + addedSteps;
    
    footprint_storage3D.clear();
    // Allocate the updated footprint vector
    footprint_storage3D.assign(desiredNumRange, vector<vector<vector<Point> > >());
    for (int t = 0; t < desiredNumRange; t++) {
        footprint_storage3D[t].assign(num_bins,vector<vector<Point> >()); }
    
    // Copy the old info into the new vector at the appropriate time step
    for (int t = 0; t < oldNumRange; t++) {
        footprint_storage3D[addedSteps + t] = oldFootprint_storage3D[t]; } 
    
    // Put the old one into the new one at the right place
    
    // Update relevant parameters
    footprint_num_range = desiredNumRange;
    // deltaT and initialUTC should remain the same
    
    
}


void Footprint3D::ShiftFootprintByMinutes(int shiftHowManyMinutes){
    // Create a temp footprint
    Footprint3D shiftedFootprint;
    shiftedFootprint = *this;    //hope to god this makes a copy!!!
    
    
    for (int minCount = 0; minCount < shiftHowManyMinutes; minCount++){
        shiftedFootprint.AddToFootprintUTC(0, 1);
//        MergeFootprintVectors(shiftedFootprint);
    }
//        this->MergeFootprintVectors(shiftedFootprint); }
//        CapeFootprint.MergeFootprintVectors(shiftedFootprint); }
    
    // Combine them all into something smooth
    SmoothedOut();
//    CapeFootprint.SmoothedOut();
    
    
    return;
}

void Footprint3D::SetFootprintUTC(double curUTC) {
    footprint_UTC = curUTC;
    UTC_Initial = curUTC;

    return;
}

double Footprint3D::GetFootprintUTC(){
    return footprint_UTC;
}





// This function assumes a constant timestep delta_t
void Footprint3D::exportGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min){
    KmlFactory* factory(KmlFactory::GetFactory());
    kmldom::DocumentPtr document(factory->CreateDocument());
    document->set_description("This is TJC description right here");
    
    
    
    // --------------------------------------------------------
    //Get all the kml stuff set up ----------------------------------------
//    KmlFactory* factory(KmlFactory::GetFactory());
//    kmldom::DocumentPtr document(factory->CreateDocument());
//    document->set_description("This is TJC description right here");
    
    
    string rocketStyleString = "m_ylw-pushpin";
    string rocketString = "s_ylw-pushpin";
    string rocketString_hl = "s_ylw-pushpin_hl";
    string rocketHref = "http://www.clker.com/cliparts/5/a/8/7/12375609571200265874pitr_Rocket_icon.svg.med.png";
    
    kmlbase::Color32 lineStyleColor;
    lineStyleColor.set_color_abgr("ff1e17ff");
    
    kmlbase::Color32 polyStyleColor;
    polyStyleColor.set_color_abgr("ff2427ff");
    
    // Create the StyleMaps -----------------------------------------||
    
    // Rocket Normal Pair     --------------------
    kmldom::PairPtr rocketPair( factory->CreatePair());
    rocketPair->set_key(0);
    rocketPair->set_styleurl(rocketString);
    
    // Rocket Highlighted Pair --------------------
    kmldom::PairPtr rocketPair_hl( factory->CreatePair());
    rocketPair_hl->set_key(1);
    rocketPair_hl->set_styleurl(rocketString_hl);
    
    kmldom::StyleMapPtr rocketStyleMap(factory->CreateStyleMap());
    rocketStyleMap->set_id(rocketStyleString);
    rocketStyleMap->add_pair(rocketPair);
    rocketStyleMap->add_pair(rocketPair_hl);
    
    document->add_styleselector(rocketStyleMap);
    
    
    
    
    
    
    // Define the style created in the StyleMap for ROCKET -------------------
//    kmldom::IconStyleIconPtr rocketIcon( factory->CreateIconStyleIcon() );
//    rocketIcon->set_href(rocketHref);
//    
//    kmldom::IconStylePtr rocketIconStylePtr( factory->CreateIconStyle() );
//    rocketIconStylePtr->set_icon(rocketIcon);
    
    kmldom::LineStylePtr lineStylePtr ( factory->CreateLineStyle() );
    lineStylePtr->set_color(lineStyleColor);
    
    kmldom::PolyStylePtr polyStylePtr ( factory->CreatePolyStyle() );
    polyStylePtr->set_color(polyStyleColor);
    
    
    kmldom::StylePtr rocketStyle( factory->CreateStyle() );
    rocketStyle->set_id(rocketString);
    rocketStyle->set_linestyle(lineStylePtr);
    rocketStyle->set_polystyle(polyStylePtr);
    //    rocketStyle->set_iconstyle(rocketIconStylePtr);

    document->add_styleselector(rocketStyle);
    
    
    // Define the style created in the StyleMap for ROCKET_HL -------------------
//    kmldom::IconStyleIconPtr rocketIcon_hl( factory->CreateIconStyleIcon() );
//    rocketIcon_hl->set_href(rocketHref);
//    
//    kmldom::IconStylePtr rocketIconStylePtr_hl( factory->CreateIconStyle() );
//    rocketIconStylePtr_hl->set_icon(rocketIcon_hl);
//    
//    kmldom::StylePtr rocketStyle_hl( factory->CreateStyle() );
//    rocketStyle_hl->set_id(rocketString_hl);
//    rocketStyle_hl->set_iconstyle(rocketIconStylePtr_hl);
//    
//    document->add_styleselector(rocketStyle_hl);
    
    

    
    // -------------------------------------------------------
    
    
    
    
    // Create the initial time structure
    time_t rawtime;
    time(&rawtime);
    
    struct tm * timeinfo;
    timeinfo = localtime ( &rawtime );  //Have to initialize the structure to the current local time then change it
    
    yyyy = yyyy - 1900; //put year into proper format
    mm = mm - 1;        //how many months SINCE january (so jan mm = 0)
    
    timeinfo->tm_year = yyyy;
    timeinfo->tm_mon = mm;
    timeinfo->tm_mday = dd;
    timeinfo->tm_hour = hour;
    timeinfo->tm_min = min;
    timeinfo->tm_sec = 0;
    
    char buffer [80];
    strftime (buffer,80, "%FT%XZ", timeinfo);
//    cout << buffer << endl;
    
    //turn rawtime into the running local time
    rawtime = mktime(timeinfo);
    
    char startTimeBuf[80];
    char endTimeBuf[80];
    
    /////
    double KmToMeters = 1e3;
	
	double alt_base = 0;
	int counter = 0;
	// t = timestep, z = z_bin, h = hull, i = a point
	for (int t = 0; t < footprint_num_range; t++) {
        
//        cout << "footprint_num_range = " << footprint_num_range << endl;
//        cout << "exporting t = " << t << endl;
        
        //Update runningTime as necessary
        time_t startTime = rawtime + t*footprint_delta_t;
        timeinfo = localtime ( &startTime );
        strftime (startTimeBuf,80, "%FT%XZ", timeinfo);
        
        time_t stopTime = rawtime + (t+1)*footprint_delta_t;
        timeinfo = localtime ( &stopTime );
        strftime (endTimeBuf,80, "%FT%XZ", timeinfo);
        
        kmldom::TimeSpanPtr timespan( factory->CreateTimeSpan() );
        timespan->set_begin(startTimeBuf);
        timespan->set_end(endTimeBuf);
        
        // Allocate the shape container
        kmldom::MultiGeometryPtr multiGeometry(factory->CreateMultiGeometry());
		
		int curr_num_bins = (INTxx) footprint_storage3D[t].size();
		for (int z=0; z < curr_num_bins; z++){
			int cur_alt_lo = (INTxx) floor((alt_base + bin_size*z)*KmToMeters);
			int cur_alt_hi = (INTxx) floor((alt_base + bin_size*(z+1))*KmToMeters);
			
			int curr_num_hulls = (INTxx) footprint_storage3D[t][z].size();
			for (int h=0; h < curr_num_hulls; h++){
				
				int curr_num_points = (INTxx) footprint_storage3D[t][z][h].size();
				for (int ix = 0; ix < curr_num_points; ix++) {
                    CoordinatesPtr coordinates(factory->CreateCoordinates());
                    
                    int firstPt = ix;
                    int secondPt = ix+1;
                    if (ix == (curr_num_points - 1)){ secondPt = 0; }
                    
                    //                    double xm1 = footprint_storage3D[t][z][h][firstPt].get_x();
                    //					double ym1 = footprint_storage3D[t][z][h][firstPt].get_y();
                    //
                    //                    double xm2 = footprint_storage3D[t][z][h][secondPt].get_x();
                    //					double ym2 = footprint_storage3D[t][z][h][secondPt].get_y();
                    //
                    //					//when using oblateness, will need to store the local_R in order to find the correct gdlat values here
                    //					double Lat1 = (180*ym1)/(PI*footprint_storage3D[t][z][h][firstPt].get_R_local());
                    //					double Lon1 = (180*xm1)/(PI*R_equator);
                    //
                    //					double Lat2 = (180*ym2)/(PI*footprint_storage3D[t][z][h][secondPt].get_R_local());
                    //					double Lon2 = (180*xm2)/(PI*R_equator);
                    
                    double Lat1 = footprint_storage3D[t][z][h][firstPt].get_gdLatDeg();
                    double Lon1 = footprint_storage3D[t][z][h][firstPt].get_lonDeg();
                    
                    double Lat2 = footprint_storage3D[t][z][h][secondPt].get_gdLatDeg();
                    double Lon2 = footprint_storage3D[t][z][h][secondPt].get_lonDeg();
                    
                    
                    coordinates->add_latlngalt(Lat1, Lon1, cur_alt_lo);       //first point, base alt
                    coordinates->add_latlngalt(Lat2, Lon2, cur_alt_lo);      //second point, base alt
                    coordinates->add_latlngalt(Lat2, Lon2, cur_alt_hi);     //second point, top alt
                    coordinates->add_latlngalt(Lat1, Lon1, cur_alt_hi);    //first point, top alt
                    
                    kmldom::LinearRingPtr linearRing(factory->CreateLinearRing());
                    linearRing->set_coordinates(coordinates);
                    kmldom::OuterBoundaryIsPtr outerBound(factory->CreateOuterBoundaryIs());
                    outerBound->set_linearring(linearRing);
                    
                    // <Polygon>
                    kmldom::PolygonPtr polygon(factory->CreatePolygon());
                    polygon->set_outerboundaryis(outerBound);
                    polygon->set_altitudemode(kmldom::ALTITUDEMODE_RELATIVETOGROUND);
                    
                    multiGeometry->add_geometry(polygon);
				}
                
                // Make the top and bottom  //Probably not necessary, could save some space here.
                // Actually just making the TOP
                {   // Using brackets so these things fall out of scope
                    CoordinatesPtr coordinates(factory->CreateCoordinates());
                    
                    for (int jx = 0; jx < curr_num_points; jx++) {
                        //                        double xm1 = footprint_storage3D[t][z][h][jx].get_x();
                        //                        double ym1 = footprint_storage3D[t][z][h][jx].get_y();
                        //                        double Lat1 = (180*ym1)/(PI*footprint_storage3D[t][z][h][jx].get_R_local());
                        //                        double Lon1 = (180*xm1)/(PI*R_equator);
                        
                        double Lat1 = footprint_storage3D[t][z][h][jx].get_gdLatDeg();
                        double Lon1 = footprint_storage3D[t][z][h][jx].get_lonDeg();
                        
                        coordinates->add_latlngalt(Lat1, Lon1, cur_alt_hi); }
                    
                    kmldom::LinearRingPtr linearRing(factory->CreateLinearRing());
                    linearRing->set_coordinates(coordinates);
                    kmldom::OuterBoundaryIsPtr outerBound(factory->CreateOuterBoundaryIs());
                    outerBound->set_linearring(linearRing);
                    
                    // <Polygon>
                    kmldom::PolygonPtr polygon(factory->CreatePolygon());
                    polygon->set_outerboundaryis(outerBound);
                    polygon->set_altitudemode(kmldom::ALTITUDEMODE_RELATIVETOGROUND);
                    
                    multiGeometry->add_geometry(polygon);
                }
                
                
				counter++;
			}
			
		}
        
        // I Think this is NOT a memory leak because, even though allocating new memory every time, I'm not losing
        //track of where it is, I'm simply storing that pointer in the folder a few lines later.
        PlacemarkPtr placemarkTemp(factory->CreatePlacemark());
        placemarkTemp->set_geometry(multiGeometry);
        placemarkTemp->set_timeprimitive(timespan);
        placemarkTemp->set_styleurl("#" + rocketStyleString);

        
        document->add_feature(placemarkTemp);
        
	}
    
    // Create <kml> and give it <Placemark>.
    kmldom::KmlPtr kml = factory->CreateKml();
    kml->set_feature(document);
    
    // Serialize to XML
    //    std::string xml = kmldom::SerializePretty(kml);
    
    // Give the xml_file a header (not sure how else to do this at the moment)
    std::string errors;
    std::string xml_output;
    kmlengine::KmlFilePtr kml_file = kmlengine::KmlFile::CreateFromImport(kml);
    kml_file->SerializeToString(&xml_output);
    
    // Look at the output in the terminal if you want
    // cout << xml_output << endl;
    
    // Write the xml to file
    bool status = kmlbase::File::WriteStringToFile(xml_output, googleEarthFile);
//    cout << "status = " << status << endl;
    
    
    
    
    return;
}


void Footprint3D::exportGoogleEarth(char *googleEarthFile){
    
    exportGoogleEarth(googleEarthFile, 2013, 5, 27, 12, 35);

    return;
}





















// This function assumes a constant timestep delta_t
//void Footprint3D::exportGoogleEarth(char *googleEarthFile){
////    static const unsigned int kHowManyPoints = 3;
//    KmlFactory* factory(KmlFactory::GetFactory());
//    
//    // Create the initial time structure
//    time_t rawtime;
//    time(&rawtime);
//    
//    struct tm * timeinfo;
//    timeinfo = localtime ( &rawtime );  //Have to initialize the structure to the current local time then change it
//    
//    int yyyy = 2013 - 1900;
//    int mm = 4; //May = 4 months since january
//    int dd = 27; //day of month
//    
//    timeinfo->tm_year = yyyy;
//    timeinfo->tm_mon = mm;
//    timeinfo->tm_mday = dd;
//    timeinfo->tm_hour = 12;
//    timeinfo->tm_min = 35;
//    timeinfo->tm_sec = 0;
//    
//    char buffer [80];
//    strftime (buffer,80, "%FT%XZ", timeinfo);
//    cout << buffer << endl;
//    
//    //turn rawtime into the running local time
//    rawtime = mktime(timeinfo);
//    cout << "raw runnign time = " << rawtime << endl;
//    
//    kmldom::DocumentPtr document(factory->CreateDocument());
//    document->set_description("This is TJC description right here");
//    
//    char startTimeBuf[80];
//    char endTimeBuf[80];
//    
//    /////
//    double KmToMeters = 1e3;
//	
//	double alt_base = 0;
//	int counter = 0;
//	// t = timestep, z = z_bin, h = hull, i = a point
//	for (int t = 0; t < footprint_num_range; t++) {
//        
//        //Update runningTime as necessary
//        time_t startTime = rawtime + t*footprint_delta_t*60.;
//        timeinfo = localtime ( &startTime );
//        strftime (startTimeBuf,80, "%FT%XZ", timeinfo);
//        
//        time_t stopTime = rawtime + (t+1)*footprint_delta_t*60.;
//        timeinfo = localtime ( &stopTime );
//        strftime (endTimeBuf,80, "%FT%XZ", timeinfo);
//        
//        kmldom::TimeSpanPtr timespan( factory->CreateTimeSpan() );
//        timespan->set_begin(startTimeBuf);
//        timespan->set_end(endTimeBuf);
//        
//        // Allocate the shape container
//        kmldom::MultiGeometryPtr multiGeometry(factory->CreateMultiGeometry());
//        
//        //let's just assume for the moment that every timestep is 10 minutes long
//        
////		int start_time = base_time + t*delta_time;		// in minutes
////		if (t == 0) { start_time = base_time - offsetTimeMinutes;}			//update start_time if on first SUA
////		int stop_time = base_time + (t+1)*delta_time;
//		
//		int curr_num_bins = (INTxx) footprint_storage3D[t].size();
//		for (int z=0; z < curr_num_bins; z++){
//			int cur_alt_lo = (INTxx) floor((alt_base + bin_size*z)*KmToMeters);
//			int cur_alt_hi = (INTxx) floor((alt_base + bin_size*(z+1))*KmToMeters);
//			
//			int curr_num_hulls = (INTxx) footprint_storage3D[t][z].size();
//			for (int h=0; h < curr_num_hulls; h++){
//				
//				int curr_num_points = (INTxx) footprint_storage3D[t][z][h].size();
//				for (int ix = 0; ix < curr_num_points; ix++) {
//                    CoordinatesPtr coordinates(factory->CreateCoordinates());
//                    
//                    int firstPt = ix;
//                    int secondPt = ix+1;
//                    if (ix == (curr_num_points - 1)){ secondPt = 0; }
//                    
////                    double xm1 = footprint_storage3D[t][z][h][firstPt].get_x();
////					double ym1 = footprint_storage3D[t][z][h][firstPt].get_y();
////                    
////                    double xm2 = footprint_storage3D[t][z][h][secondPt].get_x();
////					double ym2 = footprint_storage3D[t][z][h][secondPt].get_y();
////					
////					//when using oblateness, will need to store the local_R in order to find the correct gdlat values here
////					double Lat1 = (180*ym1)/(PI*footprint_storage3D[t][z][h][firstPt].get_R_local());
////					double Lon1 = (180*xm1)/(PI*R_equator);
////                    
////					double Lat2 = (180*ym2)/(PI*footprint_storage3D[t][z][h][secondPt].get_R_local());
////					double Lon2 = (180*xm2)/(PI*R_equator);
//                    
//                    double Lat1 = footprint_storage3D[t][z][h][firstPt].get_gdLatDeg();
//                    double Lon1 = footprint_storage3D[t][z][h][firstPt].get_lonDeg();
//                    
//                    double Lat2 = footprint_storage3D[t][z][h][secondPt].get_gdLatDeg();
//                    double Lon2 = footprint_storage3D[t][z][h][secondPt].get_lonDeg();
//                    
//                    
//                    coordinates->add_latlngalt(Lat1, Lon1, cur_alt_lo);       //first point, base alt
//                    coordinates->add_latlngalt(Lat2, Lon2, cur_alt_lo);      //second point, base alt
//                    coordinates->add_latlngalt(Lat2, Lon2, cur_alt_hi);     //second point, top alt
//                    coordinates->add_latlngalt(Lat1, Lon1, cur_alt_hi);    //first point, top alt
//                    
//                    kmldom::LinearRingPtr linearRing(factory->CreateLinearRing());
//                    linearRing->set_coordinates(coordinates);
//                    kmldom::OuterBoundaryIsPtr outerBound(factory->CreateOuterBoundaryIs());
//                    outerBound->set_linearring(linearRing);
//                    
//                    // <Polygon>
//                    kmldom::PolygonPtr polygon(factory->CreatePolygon());
//                    polygon->set_outerboundaryis(outerBound);
//                    polygon->set_altitudemode(kmldom::ALTITUDEMODE_RELATIVETOGROUND);
//                    
//                    multiGeometry->add_geometry(polygon);
//				}
//                
//                // Make the top and bottom  //Probably not necessary, could save some space here.
//                // Actually just making the TOP
//                {   // Using brackets so these things fall out of scope
//                    CoordinatesPtr coordinates(factory->CreateCoordinates());
//                    
//                    for (int jx = 0; jx < curr_num_points; jx++) {
////                        double xm1 = footprint_storage3D[t][z][h][jx].get_x();
////                        double ym1 = footprint_storage3D[t][z][h][jx].get_y();
////                        double Lat1 = (180*ym1)/(PI*footprint_storage3D[t][z][h][jx].get_R_local());
////                        double Lon1 = (180*xm1)/(PI*R_equator);
//                        
//                        double Lat1 = footprint_storage3D[t][z][h][jx].get_gdLatDeg();
//                        double Lon1 = footprint_storage3D[t][z][h][jx].get_lonDeg();
//                        
//                        coordinates->add_latlngalt(Lat1, Lon1, cur_alt_hi); }
//                    
//                    kmldom::LinearRingPtr linearRing(factory->CreateLinearRing());
//                    linearRing->set_coordinates(coordinates);
//                    kmldom::OuterBoundaryIsPtr outerBound(factory->CreateOuterBoundaryIs());
//                    outerBound->set_linearring(linearRing);
//                    
//                    // <Polygon>
//                    kmldom::PolygonPtr polygon(factory->CreatePolygon());
//                    polygon->set_outerboundaryis(outerBound);
//                    polygon->set_altitudemode(kmldom::ALTITUDEMODE_RELATIVETOGROUND);
//                    
//                    multiGeometry->add_geometry(polygon);
//                }
//                
//                
//                
//                
//				counter++;
//			}
//			
//		}
//        
//        // I Think this is NOT a memory leak because, even though allocating new memory every time, I'm not losing
//        //track of where it is, I'm simply storing that pointer in the folder a few lines later.
//        PlacemarkPtr placemarkTemp(factory->CreatePlacemark());
//        placemarkTemp->set_geometry(multiGeometry);
//        placemarkTemp->set_timeprimitive(timespan);
//        
//        document->add_feature(placemarkTemp);
//
//	}
//    
//    // Create <kml> and give it <Placemark>.
//    kmldom::KmlPtr kml = factory->CreateKml();
//    kml->set_feature(document);
//    
//    // Serialize to XML
//    //    std::string xml = kmldom::SerializePretty(kml);
//    
//    // Give the xml_file a header (not sure how else to do this at the moment)
//    std::string errors;
//    std::string xml_output;
//    kmlengine::KmlFilePtr kml_file = kmlengine::KmlFile::CreateFromImport(kml);
//    kml_file->SerializeToString(&xml_output);
//    
//    // Look at the output in the terminal if you want
//    // cout << xml_output << endl;
//    
//    // Write the xml to file
//    bool status = kmlbase::File::WriteStringToFile(xml_output, googleEarthFile);
//    cout << "status = " << status << endl;
//    
//    
//    
//    
//    return;
//}





// ~~~~~~~~~~~~~~~~~~~~~~~~~ !!!!GRAVEYARD!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//void Footprint3D::create_a_footprint(int cur_index) {
//	
//	// bin the points based on their z-location
//	int num_current_points = all_points.size();
//    
//	for (int i = 0; i < num_current_points; i++) {
//		int bin_index = floor(all_points[i].get_z()/bin_size);
//		if ((bin_index >=0) && (bin_index < num_bins)) {
//			current_points[bin_index].push_back(all_points[i]); }			//BAD ACCESS ERROR HERE
//        //Think it comes from hardcoding the min and max bin values
//		else {
//			// If you wound up here, that means the altitude of the current point is out of the range you specified for the bins.
//		} }
//	
//	cout << "done binning by altitude" << endl;
//    
//	
//    //	cout << "Let's see if I did that correctly!\n";
//    //	int test_index = 5;
//    //	for (int i = 0; i < current_points[test_index].size(); i++) {
//    //		cout << current_points[test_index][i] << "   " << current_points[test_index][i].get_id() << endl; }
//	
//	
//    
//	//problem bin is 16315
//	for (int cur_bin = 0; cur_bin < num_bins; cur_bin++) {
//		//sort the points, keeping only the unique ones
//        //	vector<Point> current_points = all_points[cur_index];
//		
//        //		cout << "start sorting.  current_points[bin] MB = " << current_points[cur_bin].capacity()*sizeof(Point)*bytes2mega << endl;
//        
//        //		cout << cur_bin << ", ";
//		
//		vector<Point> unique_points = sort_unique_points(cur_bin);
//		int num_unique_pts = unique_points.size();
//		
//        //		cout << "unique points:\n";
//        //		for (int jx = 0; jx < num_unique_pts; jx++){
//        //			cout << unique_points[jx] << endl; }
//		
//        //		cout << "done sorting.  unique_points MB = " << unique_points.capacity()*sizeof(Point)*bytes2mega << endl;
//		
//		//----------------------- Swinging Arm -------------------------
//		
//		
//		
//		
//		vector<int> hull_indices;
//		
//		int vec_size = num_unique_pts;
//		
//		vector<bool> available;
//		available.assign(vec_size,true);	//make all pts available
//		
//		vector<int> avail_IX;				//available indices
//		
//		vector<bool> visited;
//		visited.assign(vec_size, false);	//set all pts as unvisited
//		
//		int num_available;
//		
//		vector<Point> segments;
//		vector<Point> temp_segments;
//        
//		for (int i = 0; i < vec_size; i++) {
//			if (available[i]) {
//				avail_IX.push_back(i); } }
//		num_available = avail_IX.size();
//		
//		int cur_hull = 0;
//		while (num_available != 0) {
//			hull_indices.clear();	//gets reset for every sub-hull
//			hull_indices.push_back(avail_IX[0]);	//places highest y-point as first hull point (index)
//			
//			double angle_offset = 0;
//			double last_angle = 0;
//			
//			int j=0;
//			bool alive = true;
//			
//			
//			int cur_pt;
//			
//			while (alive && ( j < 2*num_available ) ) {
//				cur_pt = hull_indices.back();		//gets the most recent hull point (index)
//				
//				double last_best_angle = 0;
//				double last_best_pt = cur_pt;
//				
//				//check to see which points are available
//				//this kind of sucks, pull out the j==0 condition to before the loop, and do the else at the end of loop
//				//	and make it a function
//				
//				
//				double x1 = unique_points[cur_pt].get_x();
//				double y1 = unique_points[cur_pt].get_y();
//				
//				//look through all available points for the one that's within arms length and
//				//  forms the smallest angle from the previous line segment
//				for (int k = 0; k < num_available; k++) {
//					double x2 = unique_points[avail_IX[k]].get_x();
//					double y2 = unique_points[avail_IX[k]].get_y();
//					
//					double dx = x2-x1;
//					double dy = y2-y1;
//					
//					double range = sqrt(dx*dx + dy*dy);
//					
//					if (cur_pt == avail_IX[k]) {	//throw out current point from calculations
//						range = 1e6; }
//					
//					if ((range > 0.0) && (range < arm)) {
//						//get angle and transform it be be from [0,2pi] measured clockwise from +x axis
//						double angle = atan2(dy, dx);
//						if (dy < 0) {
//							angle = -angle; }
//						else {
//							angle = 2*PI - angle; }
//						
//						
//						//some angles will appear to be smaller (better), but really they're nearly 360 degrees
//						//	away from the last segment.  This catches those guys, as well as 360 deg rotations.
//						//NOTE: MAY BE A BUG HERE!!! Only catches 360 rots on the upswing (offset neg) ???
//						if (((fabs(angle)-eps) <= fabs(angle_offset)) && (angle_offset < 0)) {
//							angle = angle + 2*PI; }
//						
//						double angle_from_previous_line = angle + angle_offset;
//						if (((angle_from_previous_line) < last_best_angle) || (last_best_angle == 0)) {
//							
//							//add in the current segment
//							temp_segments.push_back(unique_points[cur_pt]);
//							temp_segments.push_back(unique_points[avail_IX[k]]);
//							
//							//do they intersect?
//							int num_intersections = do_they_intersect(&temp_segments);
//							
//							//remove the current segment so we can check for others on the next iteration
//							temp_segments.pop_back();
//							temp_segments.pop_back();
//							
//							if ((num_intersections == 0) || ((num_intersections == 1) && (avail_IX[k] == hull_indices[0]))){
//								last_best_angle = angle_from_previous_line;
//								last_best_pt = avail_IX[k];
//								last_angle = angle; }
//							else { /*this isn't allowed, ignore it*/ } }
//						
//						
//					} //ends range if
//					
//					
//				} //ends k loop
//				
//				//save the good stuff
//				segments.push_back(unique_points[hull_indices.back()]);
//				segments.push_back(unique_points[last_best_pt]);
//				hull_indices.push_back(last_best_pt);
//				
//				//this vector isnt actually necessary, but i think it makes me a little safer
//				temp_segments.clear();
//				temp_segments = segments;
//				
//				visited[hull_indices[j+1]] = true;
//				
//				angle_offset = PI - last_angle;
//				
//				if (hull_indices[0] == hull_indices.back()) {
//					alive = false; }
//				
//				
//				
//				j++;
//			}//ends while loop that creates a single footprint
//			
//			
//			segments.clear();
//			
//            //			if (cur_bin == 4) {
//            //				cout << "entering code for debug" << endl; }
//			
//			update_available_points(&available, &avail_IX, unique_points, hull_indices);
//			num_available = avail_IX.size();
//			
//			//--------------------- Store the footprint -------------
//			
//			//NEED SOME MAGIC HERE IN ALLOCATING THE STORAGE VECTOR.  THIS WILL BE AN ERROR.
//			
//            //			footprint_storage[cur_index].push_back(vector<Point>());
//            //			for (int i = 0; i < hull_indices.size(); i++) {
//            //				footprint_storage[cur_index][cur_hull].push_back(unique_points[hull_indices[i]]); }
//			
//			footprint_storage3D[cur_index][cur_bin].push_back(vector<Point>());
//			for (int i = 0; i < hull_indices.size(); i++) {
//				footprint_storage3D[cur_index][cur_bin][cur_hull].push_back(unique_points[hull_indices[i]]); }
//			
//			
//			
//			// Start debugging statements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			
//            //			if (cur_hull == 0) {
//            //				double x_curr = current_points[cur_bin][0].get_x();
//            //				double y_curr = current_points[cur_bin][0].get_y();
//            //
//            //				cout << "all_bin_pts = [\n";
//            //
//            //				cout << std::setprecision(12) << x_curr << "  " << y_curr;
//            //
//            //				for (int i = 1; i < current_points[cur_bin].size(); i++) {
//            //					x_curr = current_points[cur_bin][i].get_x();
//            //					y_curr = current_points[cur_bin][i].get_y();
//            //
//            //					cout << ";\n" << x_curr << "  " << y_curr;
//            //				}
//            //				cout << "];\n\n";
//            //			}
//			
//            //			cout << "cur_bin = " << cur_bin << endl;
//            //
//            //			// Print the hull boundary pts
//            //			double x_hull = footprint_storage3D[cur_index][cur_bin][cur_hull][0].get_x();
//            //			double y_hull = footprint_storage3D[cur_index][cur_bin][cur_hull][0].get_y();
//            //
//            //			if (cur_hull == 0) {
//            //				cout << "curr_hull = [\n"; }
//            //			else {
//            //				cout << "curr_subhull = [\n"; }
//            //
//            //			cout << std::setprecision(12) << x_hull << "  " << y_hull;
//            //
//            //			for (int i = 1; i < hull_indices.size(); i++) {
//            //				x_hull = footprint_storage3D[cur_index][cur_bin][cur_hull][i].get_x();
//            //				y_hull = footprint_storage3D[cur_index][cur_bin][cur_hull][i].get_y();
//            //
//            //				cout << ";\n" << x_hull << "  " << y_hull;
//            //			}
//            //			cout << "];\n\n";
//            //
//            //			// Print the available points
//            //			double x_avail = unique_points[avail_IX[0]].get_x();
//            //			double y_avail = unique_points[avail_IX[0]].get_y();
//            //
//            //			if (avail_IX.size() > 0) {
//            //				cout << "avail" << cur_hull << " = [\n";
//            //				cout << std::setprecision(12) << x_avail << "  " << y_avail;
//            //
//            //				for (int i = 1; i < avail_IX.size(); i++) {
//            //					x_avail = unique_points[avail_IX[i]].get_x();
//            //					y_avail = unique_points[avail_IX[i]].get_y();
//            //					
//            //					cout << ";\n" << x_avail << "  " << y_avail;
//            //				}
//            //				cout << "];\n\n";
//            //				
//            //				
//            //				for (int i = 0; i < avail_IX.size(); i++) {
//            //					cout << avail_IX[i] << "  "; }
//            //				cout << endl << endl;
//            //			}
//			
//			
//			// End debugging statements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			
//			
//			
//			
//			
//			cur_hull++;
//			
//            //			for (int i = 0; i < hull_indices.size(); i++) {
//            //				boundary_point_storage[cur_index].push_back(unique_points[hull_indices[i]]); }
//			
//			
//			
//		}//ends outer while loop that generates all footprints
//        
//        
//        //--------------------- End Swinging Arm -------------------------
//        
//        
//        unique_points.clear();
//		
//	} //ends the loop over the bins
//}



//void Footprint3D::MergeFootprintVectors(string new_footprint_name, string FPOut_name){
//
//	// -------- Copy the vitals of the existing footprint (call it FP2)
//	int FP2_footprint_num_range = footprint_num_range;
//	double FP2_footprint_UTC = footprint_UTC;
//	double FP2_footprint_delta_t = footprint_delta_t;
//	double FP2_bin_size = bin_size;
//	int FP2_num_bins = num_bins;
//
////	// Not vital, but will want to use it.
////	double FP2_UTC_Final = UTC_Final;
//	double FP2_finalUTC = FP2_footprint_UTC + FP2_footprint_num_range*FP2_footprint_delta_t/(24.*60.);
//
//
//
//	vector<vector<vector<vector<Point> > > > FP2_storage3D;	//[cur_index][z_bin][cur_hull][point]
//	FP2_storage3D = footprint_storage3D;
//
//	// --------- Wipe out the old footprint and load up the new one
//	load_footprint_as_vector(new_footprint_name);
//
//	// Calculate a few things that the main constructor calculates just to be consistent
//	double UTC_Initial = footprint_UTC;
//	double delta_t = footprint_delta_t;
//	double UTC_Final = UTC_Initial + delta_t*footprint_num_range/(60*24);
//
//	// What's the total range?
//	int TEMP_num_range = (INTxx) round(fabs(UTC_Final - FP2_finalUTC)*24*60/delta_t);
//	cout << "tempnumrange = " << TEMP_num_range << endl;
//
//
//	if (UTC_Final > FP2_finalUTC) {
//		TEMP_num_range = footprint_num_range; }
//	else {
//		TEMP_num_range = FP2_footprint_num_range; }
//
//	// ---------- Now start to combine them into a single temp_vector
//	int FP1start = -2;
//	int FP2start = -2;
//
//	if (footprint_UTC == FP2_footprint_UTC){
//		cout << "The UTC's are equal!" << endl; }
//	else if (footprint_UTC > FP2_footprint_UTC) {
//		// This means that FP2 starts first
//		FP2start = 0;
//
//		// How many delta_t's before FP1 starts?
//		cout << "fp1start = " << (footprint_UTC - FP2_footprint_UTC)*24*60 / delta_t << endl;
//		FP1start = (INTxx) round((footprint_UTC - FP2_footprint_UTC)*24*60 / delta_t);
//		cout << "again but as an int = " << FP1start << endl;
//
//	}
//	else if (footprint_UTC < FP2_footprint_UTC) {
//		// This means that FP1 starts first
//		FP1start = 0;
//
//		// How many delta_t's before FP2 starts?
//		cout << "fp2start = " << (FP2_footprint_UTC - footprint_UTC)*24*60 / delta_t << endl;
//		FP2start = (INTxx) round((FP2_footprint_UTC - footprint_UTC)*24*60 / delta_t);
//		cout << "again but as an int = " << FP2start << endl;
//	}
//	else {
//		cout << "The UTC's for the two footprints aren't equal, this case isn't handled yet";
//		exit(9); }
//
//	double TEMP_num_bins = -5;
//	if (num_bins == FP2_num_bins){
//		cout << "The number of bins are equal!" << endl;
//		TEMP_num_bins = num_bins; }
//	else {
//		cout << "The bins for the two footprints aren't equal, this case isn't handled yet";
//		exit(9); }
//
//
//
//	vector<vector<vector<vector<Point> > > > TEMP_storage3D;	//[cur_index][z_bin][cur_hull][point]
//	TEMP_storage3D.assign(TEMP_num_range, vector<vector<vector<Point> > >());
//
//
//
//	// Do the merge
//	for (int t = 0; t < TEMP_num_range; t++) {
//		TEMP_storage3D[t].assign(TEMP_num_bins, vector<vector<Point> >());
//
//		// These indices will keep track of the separate footprints
//		int t1 = t - FP1start;
//		int t2 = t - FP2start;
//
//		for (int binx = 0; binx < FP2_num_bins; binx++) {
//			int existing_hulls = 0;
//			int merging_hulls = 0;
//
//			if ((t1 >= 0) && (t1 < footprint_num_range)) {
//				existing_hulls = (INTxx) footprint_storage3D[t1][binx].size(); }
//
//			if ((t2 >= 0) && (t2 < FP2_footprint_num_range)){
//				merging_hulls = (INTxx) FP2_storage3D[t2][binx].size(); }
//
//			int num_hulls = merging_hulls + existing_hulls;
//			TEMP_storage3D[t][binx].assign(num_hulls, vector<Point>());
//
//			int curr_hull = 0;
//			for (int hullx = 0; hullx < existing_hulls; hullx++) {
//				int num_points = (INTxx) footprint_storage3D[t1][binx][hullx].size();
//				cout << "num poins = " << num_points << endl;
//				TEMP_storage3D[t][binx][curr_hull].assign(num_points, Point());
//				memcpy((void *) &TEMP_storage3D[t][binx][curr_hull][0], (void *) &footprint_storage3D[t1][binx][hullx][0], num_points*sizeof(Point));
//				curr_hull++; }
//
//			for (int hullx = 0; hullx < merging_hulls; hullx++) {
//				int num_points = (INTxx) FP2_storage3D[t2][binx][hullx].size();
//				cout << "num pointtts = " << num_points << endl;
//				TEMP_storage3D[t][binx][curr_hull].assign(num_points, Point());
//				memcpy((void *) &TEMP_storage3D[t][binx][curr_hull][0], (void *) &FP2_storage3D[t2][binx][hullx][0], num_points*sizeof(Point));
//				curr_hull++; }
//		}
//	}
//
//	// ----- Wipe out existing footprint and replace with new one
//	Destruct_footprint();
//
//	// -------- Copy the vitals of the footprint
//	footprint_num_range = TEMP_num_range;
//	footprint_UTC = FP2_footprint_UTC;
//	footprint_delta_t = FP2_footprint_delta_t;
//	bin_size = FP2_bin_size;
//	num_bins = FP2_num_bins;
//
//	footprint_storage3D = TEMP_storage3D;
//
//	// Calculate a few things that the main constructor calculates just to be consistent
//	UTC_Initial = FP2_footprint_UTC;
//	delta_t = FP2_footprint_delta_t;
//	UTC_Final = UTC_Initial + delta_t*TEMP_num_range/(60*24);
//
//	// Write the merged vector to a file (if asked to)
//	store_footprint_as_vector(FPOut_name);
//
//	cout << "leaving merge" << endl;
//
//	return;
//}





//void Footprint3D::store_footprint_as_points(){
//// It's not clear to my why I would ever use this function...for appending, should always write the footprint_storage3D vector
//
//	// Get the points into the appropriate vector
//	std::vector<std::vector<Point> > total_points_at;	//trajectory and debris points all put together by timestep
//	total_points_at.assign(footprint_num_range,std::vector<Point>());
//
//	for (int t = 0; t < footprint_num_range; t++){
//		for (int ix = 0; ix < footprint_storage3D[t].size(); ix++) {
//			for (int jx = 0; jx < footprint_storage3D[t][ix].size(); jx++) {
//				for (int kx = 0; kx < footprint_storage3D[t][ix][jx].size(); kx++) {
//					total_points_at[t].push_back(footprint_storage3D[t][ix][jx][kx]); }}}}
//
//
//
//	// Write the output file used to make a footprint out of
//	ofstream outfile;
//	outfile.open("GeneratedFiles/points_to_wrap_up_abridged.dat",ios::out | ios::binary);
//
//	outfile.write((char *) &footprint_num_range, sizeof(footprint_num_range));
//	outfile.write((char *) &footprint_UTC, sizeof(footprint_UTC));
//	outfile.write((char *) &footprint_delta_t, sizeof(footprint_delta_t));
//	//	cout << "put ptr = " << outfile.tellp() << endl;
//	//	outfile.write((char *) &num_debris_events, sizeof(num_debris_events));
//
//	int num_points_here;
//	for (int t = 0; t < footprint_num_range; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here));
//	}
//
//	for (int t = 0; t < footprint_num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point));
//	}
//
//	outfile.close();
//
//	return;
//}















