/*
 *  Architecture.cpp
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 7/19/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 *
 *  Archs to look at
 *	- Oldschool Eastern Range (3 hours early)
 *	- Half of Eastern Range depending on launch azimuth (3 hours early and only a few minutes early)
 *	- Tubes of various sizes that only capture the nominal profile
 *	- SUA that encloses all debris positions given explosion at all timesteps
 *		These last two in conjunction can show us where there is diminishing returns in terms of SUA size
 *			and also can use size of debris area with response time of ATC and aircraft to determine the smallest
 *			tube radius that could be used.
 *	- Launch using nominal tube that has debris event which activates a debris SUA that planes must flee from.
 *		For different tube sizes and debris SUA areas, can see just how much you gain with this approach over the simpler
 *			tube-all-the-time or all-debris-all-the-time scenarios.
 *  - Not time evolving, simply debris pattern extremum pattern from ground to sky for explosion along nominal 
 *
 *	More Complicated
 *	- Simulate that holding pattern containment area thing from the FAA Conops paper 
 *	- Naturally want to construct something for Virgin Galactic
 *
 *  NOTE: Don't just look at the established ranges.  It's launch from new spaceports that is most available
 *		to be impacted by this kind of analysis.
 *
 *
 *  * Within these functions, IT IS NOT APPROPRIATE to truncate the points at the top of the NAS.  That happens in Footprint3D.
 *		The points that are output here should be more general
 *
 * */


#include "Architecture.h"

Architecture::Architecture() {

	return;
}




void Architecture::read_single_trajectory_from_file(char *trajectoryFileName){
	//Need to be careful when writing to a binary.  the elements of a std::vector are allocated such that they
	//  are memory contiguous, but the contiguity of a vector of a vector is NOT guaranteed.
	
	//	int num_to_write = 1;	//would be num_per_batch if we were writing a bunch
	
	ifstream infile;
	infile.open(trajectoryFileName,ios::in | ios::binary);
	
	std::vector<double> ZeroVec3Size(3,0.);
	
	int num_trajectories;
	if (infile.good()) {
		
		infile.read((char *) &num_trajectories, sizeof(num_trajectories));	// read in how many time steps there are
		RocketLatLonAltStorage.assign(num_trajectories,std::vector<std::vector<double> >() );
		cout << "num_trajectories = " << num_trajectories << endl;
		
		int num_points_here;
		for (int ix = 0; ix < num_trajectories; ix++) {
			infile.read((char *) &num_points_here, sizeof(num_points_here));
			cout << "num_point_here = " << num_points_here << endl;
			RocketLatLonAltStorage[ix].assign(num_points_here,ZeroVec3Size);}
				
		//Read in each value and assign it
		for (int ix = 0; ix < num_trajectories; ix++) {
			num_points_here = (INTxx) RocketLatLonAltStorage[ix].size();
			cout << "num_point_here = " << num_points_here << endl;
			for (int jx = 0; jx < num_points_here; jx++) {
				for (int kx = 0; kx < 3; kx++) {
					infile.read((char *) &RocketLatLonAltStorage[ix][jx][kx], sizeof(double)); } } }
		
		// Now read the basic timing information
		infile.read((char *) &Initial_UTC, sizeof(Initial_UTC));
		infile.read((char *) &num_time_steps, sizeof(num_time_steps));
		
		// Allocate the StateTimes array (would rather have this be a std::vector...try that later) //
		RocketStateTimes = new double[num_time_steps];
		
		// Read the StateTimes array (simple array of doubles, i expect this to be memory contiguous)
		infile.read((char *) &RocketStateTimes[0], num_time_steps*sizeof(double));
		
		// Read the size of the state vector
		infile.read((char *) &state_vector_size, sizeof(state_vector_size));
		std::vector<double> ZeroVecStateSize(state_vector_size,0.);

		// Now read the state vectors in ECI at every timestep
		State_Vector_Storage_Vec.assign(num_trajectories,std::vector<std::vector<double> >());
		for (int ix = 0; ix < num_trajectories; ix++) {
			num_points_here = (INTxx) RocketLatLonAltStorage[ix].size();
			State_Vector_Storage_Vec[ix].assign(num_points_here, ZeroVecStateSize);
			for (int jx = 0; jx < num_points_here; jx++) {
				for (int kx = 0; kx < state_vector_size; kx++) {
					infile.read((char *) &State_Vector_Storage_Vec[ix][jx][kx], sizeof(double)); } } } 
		
		infile.read((char *) &first_stage_struct_mass, sizeof(double));
		
	} else {
		cout << "Opening trajectory binary file failed" << endl;
	}
	
	infile.close();

	return;
}

void Architecture::read_debris_from_file(unsigned int idNum){
	string fileString = "GeneratedFiles/debris_points";
	string extString = ".dat";
	char buffer[5];
	sprintf(buffer,"%i",idNum);
	fileString = fileString + buffer + extString;
	
	ifstream infile;
	infile.open(fileString.c_str(),ios::in | ios::binary);
	
	std::vector<double> ZeroVec3Size(3,0.);
	
	int num_pieces;
	
	if (infile.good()) {
		infile.read((char *) &num_pieces, sizeof(num_pieces));	// read in how many time steps there are
		DebrisLatLonAltStorage.assign(num_pieces,std::vector<std::vector<double> >() );
		
		int num_points_here;
		for (int ix = 0; ix < num_pieces; ix++) {
			infile.read((char *) &num_points_here, sizeof(num_points_here));
			DebrisLatLonAltStorage[ix].assign(num_points_here,ZeroVec3Size);}
				
		//Read in each value and assign it
		for (int ix = 0; ix < num_pieces; ix++) {
			num_points_here = (INTxx) DebrisLatLonAltStorage[ix].size();
			for (int jx = 0; jx < num_points_here; jx++) {
				for (int kx = 0; kx < 3; kx++) {
					infile.read((char *) &DebrisLatLonAltStorage[ix][jx][kx], sizeof(double)); } } }
		
		// Now read the basic timing information
		infile.read((char *) &DebrisInitialUTC, sizeof(DebrisInitialUTC));
		infile.read((char *) &debris_num_time_steps, sizeof(debris_num_time_steps));
		
		// Allocate the StateTimes array (would rather have this be a std::vector...try that later)
		DebrisStateTimes = new double[debris_num_time_steps];
		
		// Read the StateTimes array (simple array of doubles, i expect this to be memory contiguous)
		infile.read((char *) &DebrisStateTimes[0], debris_num_time_steps*sizeof(double));
		
		// Calculate final UTC
		DebrisFinalUTC = DebrisInitialUTC + DebrisStateTimes[debris_num_time_steps-1]/(3600*24);
		
	} else {
		cout << "Opening debris file failed: " << fileString << endl;
	}
	
	infile.close();
	
	
	return;
}

void Architecture::reset_debris(){
	
	DebrisLatLonAltStorage.clear();			// THIS MAKES ME UNCOMFORTABLE
	DebrisInitialUTC = -2;
	debris_num_time_steps = -2;
	delete[] DebrisStateTimes;
	
	return;
}

double Architecture::get_local_earth_radius(double gdlat /*rad*/){
	double a = R_equator;
	double b = R_polar;
	
	return sqrt((pow(a*a*cos(gdlat),2) + pow(b*b*sin(gdlat),2))/(pow(a*cos(gdlat),2) + pow(b*sin(gdlat),2)));
}


void Architecture::make_tube_around_single_trajectory(double r_desired /*km*/, int fidelity){
	// higher fidelity means more points used to make the tube
	double delta_angle = PI/(fidelity-1);
	
	int num_points_in_trajectory = (INTxx) RocketLatLonAltStorage[0].size();
	cout << "pts in traj = " << num_points_in_trajectory << endl;
	Single_Trajectory_Tube.assign(num_points_in_trajectory, std::vector<Point>());
	Point EmptyPoint;
		
	int num_fidelity_points = 2*fidelity-2;
	for (int cur_pt = 0; cur_pt < num_points_in_trajectory; cur_pt++){
		Single_Trajectory_Tube[cur_pt].assign(num_fidelity_points,EmptyPoint);

		double gdlat0 = RocketLatLonAltStorage[0][cur_pt][0];
		double lon0 = RocketLatLonAltStorage[0][cur_pt][1];
		double alt0 = RocketLatLonAltStorage[0][cur_pt][2];
		double local_R = get_local_earth_radius(gdlat0);
		
//		cout << "gdlat = " << gdlat0*180/PI << "      local_R = " << local_R << endl;
		
		//left off here
		double x0 = R_equator*lon0;
		double y0 = local_R*gdlat0;
		
		int cur_tube_pt = 0;
		int ix = 0;

		double x = x0 + r_desired*cos(ix*delta_angle);
		double y = sqrt(r_desired*r_desired - (x-x0)*(x-x0)) + y0;
		
//		cout << "x0 = " << x0 << "   x = " << x << "  lon = " << x/R_equator * (180/PI) << endl;
//		cout << "y0 = " << y0 << "   y = " << y << "  lat = " << y/local_R * (180/PI)<< endl;
		
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_x(x);
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_y(y);
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_z(alt0);
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_R_local(local_R);
		cur_tube_pt++;
		ix++;
		
		while (ix < fidelity-1) {
			x = x0 + r_desired*cos(ix*delta_angle);
			y = sqrt(r_desired*r_desired - (x-x0)*(x-x0)) + y0;
			
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_x(x);
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_y(y);
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_z(alt0);
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_R_local(local_R);
			
			cur_tube_pt++;
			
			y = -sqrt(r_desired*r_desired - (x-x0)*(x-x0)) + y0;
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_x(x);
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_y(y);
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_z(alt0);
			Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_R_local(local_R); 
		
			cur_tube_pt++;
			ix++;
		}

		x = x0 + r_desired*cos(ix*delta_angle);
		y = sqrt(r_desired*r_desired - (x-x0)*(x-x0)) + y0;
		
		
		
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_x(x);
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_y(y);
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_z(alt0);
		Single_Trajectory_Tube[cur_pt][cur_tube_pt].set_R_local(local_R);
	}
		
	return;
}

void Architecture::explode_at_all_points() {
	char WindOption[] = "simple atmosphere";
	char DensityOption[] = "cantwell density";
	char StateOption[] = "debris";

	double debris_delta_t = 1.0;
	Debris TestCatalog;
	TestCatalog.GenerateRandomPieces();

	Trajectory MyDebris(StateOption, WindOption, DensityOption);
	
	int currTraj = 0;
	double debris_state[8];	//contains state vector at time of breakup PLUS the UTC time
	
	for (int tstep = 0; tstep < num_time_steps; tstep++) {
		// Copy the state vector into debris_state
		memcpy((void *) debris_state, (void *) &State_Vector_Storage_Vec[currTraj][tstep][0], 7*sizeof(double));
		debris_state[7] = Initial_UTC + RocketStateTimes[tstep]/(24*3600);	// And append the UTC time
		MyDebris.Propagate_Debris_From_Catalog(debris_state, TestCatalog, debris_delta_t, (unsigned int) tstep);}
	
	return;
}














vector<vector<Point> > Architecture::write_all_points_debris(double tstepMinutes, string outFileName) {
	// Writes the points assuming instantaneous.  
	// In the future, could have a lookahead time and group debris locations with that (note: that's not binning).
	
	// Dump those points into the all_points vector (initialized to accomodate all (desired) time steps
	//Collect all the points
	//Allocate the vector of points that we're going to write to file
	double timeInFlight = (DebrisStateTimes[debris_num_time_steps-1] - DebrisStateTimes[0])/60.;	//minutes
	int time_steps_out = (INTxx) ceil(timeInFlight/tstepMinutes);
	
	std::vector<std::vector<Point> > total_points_at;	
	total_points_at.assign(time_steps_out,std::vector<Point>());
	
	
//	// Open a file to write a ground signature for debugging / validation purposes
//	ofstream gsig;
//	gsig.open("/Users/marian/Documents/MATLAB/tjc_ground_debug2.m",ios::out);
//	ofstream gxyz;
//	gxyz.open("/Users/marian/Documents/MATLAB/tjc_xyz_debug.m",ios::out);
//	cout << "HEY!!!  You're secretly writing a ground signature to a file.  Delete this when done debugging~~~~~~~~~~~~~~" << endl;
//	gsig << "tjc_ground = [";
//	gxyz << "tjc_xyz = [";
	
	
	int numPieces = (INTxx) DebrisLatLonAltStorage.size();
	for (int pc = 0; pc < numPieces; pc++) {
		int numCurrTimeSteps = (INTxx) DebrisLatLonAltStorage[pc].size();
		for (int tstep = 0; tstep < numCurrTimeSteps; tstep++) {
			// Figure out the timestep index we're going to bin this into
			int binTimeStep = (INTxx) floor(DebrisStateTimes[tstep]/(60.*tstepMinutes));
			
			// I recall the lon == 0 condition fixed some problem when not all points were initialized.
			//		I don't know if this is still a problem or not.  Keep it for now.
			if (!(DebrisLatLonAltStorage[pc][tstep][1] == 0)
				&& (binTimeStep < time_steps_out)) {
				
				Point temp_pt;
				double temp_local_R = get_local_earth_radius(DebrisLatLonAltStorage[pc][tstep][0]);
				
				temp_pt.set_x(DebrisLatLonAltStorage[pc][tstep][1]*R_equator);
				temp_pt.set_y(DebrisLatLonAltStorage[pc][tstep][0]*temp_local_R);
				temp_pt.set_z(DebrisLatLonAltStorage[pc][tstep][2]);
				temp_pt.set_R_local(temp_local_R);
				
				total_points_at[binTimeStep].push_back(temp_pt); } 
		
//			if (tstep == numCurrTimeSteps-1) {
//				gsig << DebrisLatLonAltStorage[pc][tstep][0]*180/PI << "  "
//					<< DebrisLatLonAltStorage[pc][tstep][1]*180/PI << "  "
//					<< DebrisLatLonAltStorage[pc][tstep][2]*1e3  << endl;
//
//				gxyz << total_points_at[binTimeStep].back().get_x() << "  "
//				<< total_points_at[binTimeStep].back().get_y() << "  "
//				<< total_points_at[binTimeStep].back().get_z() << "  "
//					<< endl;
//			}
		
		} }
//	gsig << "];\n\n\n";
//	gsig.close();
//
//	gxyz << "];\n\n\n";
//	gxyz.close();
	
		
	if (outFileName != "none") {
		// Now write the stuff to a file
		ofstream outfile;
		outfile.open(outFileName.c_str(),ios::out | ios::binary);
		
		// Write timing info
		outfile.write((char *) &time_steps_out, sizeof(time_steps_out));
		outfile.write((char *) &DebrisInitialUTC, sizeof(DebrisInitialUTC));
		outfile.write((char *) &tstepMinutes, sizeof(tstepMinutes));
		
		// Write vector structure info
		int num_points_here;
		for (int t = 0; t < time_steps_out; t++) {
			num_points_here = (INTxx) total_points_at[t].size();
			outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
		
		// Write the points
		for (int t = 0; t < time_steps_out; t++) {
			outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
		
		outfile.close();
	}
	
	return total_points_at;
}


/*! MonteCarlo: Currently, this function
 *  * propagates a first stage trajectory
 *  * creates a very simple LRHC envelope approximation around the trajectory
 *  * creates the actual envelope around LRHC with a specified bin size 
 *  * propagates debris to the ground starting from altitude_we_explode_until_km at every tstep = tstep - backwardsStep
 *  * at each considered timestep (of initial explosion) APPENDS the generated points into the existing envelope
 *  * when finished, stores the footprint as a vector in a binary file via store_footprint_as_vector
 */
//void Architecture::MonteCarlo(char WindOption[], char DensityOption[], string nominal_traj_filename, string footprintVectorFile) {
//// NOTE: Everything will be automatically calculated to have a timestep of 1 minute.  When everything has been calculated,
////  we can then easily translate the footprint timestep into any integer multiple of 1 minute afterwards.
//	
//	// ~~~~~ Propagate the first stage trajectory ~~~~~~~
//	char FirstStage_StateOption[] = "first stage";
//	int num_launches_per_batch = 1;								// I think this might be obsolete given how the structure of the code is evolving.  Leave it at 1
//	double delta_t = 2.;										// Size of time step at which solution will be sampled (solver uses adaptive step method)
//
//	//Initialize the atmosphere and let it know to expect a first-stage
//	Trajectory MyTraj(FirstStage_StateOption, WindOption, DensityOption);
//	MyTraj.Initialize_First_Stage(num_launches_per_batch, delta_t, nominal_traj_filename);
//	MyTraj.Propagate_First_Stage();
//	MyTraj.write_single_trajectory_to_file("GeneratedFiles/trajectory_points.dat");	//Writes file to "GeneratedFiles/trajectory_points.dat"
//
//	// Propagation-related functions that aren't being used right now
//		// MyTraj.Transform_Nominal("Files/zzz_nominal_cape_finished.txt", "not currently using this field");
//		// MyTraj.write_batch_to_google_earth("useless", 0);
//
//	// ~~~~~ Generate LRHC and wrap it up ~~~~~~~~~~
//	// Want to have a starting footprint to append to, so use the LRHC envelope which we roughly simulate
//	
//	// Things needed to create LRHC	
////	string nominalFileName("GeneratedFiles/trajectory_points.dat");
//    char nominalFileName[] = "GeneratedFiles/trajectory_points.dat";
//	int fidelity = 16;
//	double semiMajorPercent = 1.5;
//	double eccentricity = 0.95;
////	SimLeftRightHotCold(nominalFileName, fidelity, semiMajorPercent, eccentricity);
//	
//	// Write the LRHC points to file
//	string lrhcFileName("GeneratedFiles/points_to_wrap_up_LRHC.dat");	//already done
//	bool isInstantaneous = true;
//	double tstepMinutes = 1;
//	double startUTC = Initial_UTC;
//	double minutesOn = 5.;
////	write_LRHC_points(lrhcFileName, isInstantaneous, tstepMinutes, startUTC, minutesOn);
//		
//	// Create the starting envelope from those LRHC points
//	double binSizeKm = 5.;
//
//	Footprint3D my_footprint(lrhcFileName,binSizeKm);
//	my_footprint.generate_footprint_at_timesteps();
//	
//    leftoffhere
//	int num_debris_runs_per_explosion = 1;
//	double altitude_we_explode_until_km = 80;
//	
//	// Set up the debris to be used
//	char Debris_StateOption[] = "debris";
//
//	double debris_delta_t = 1.0;
//	Debris TestCatalog;
//	Trajectory MyDebris(Debris_StateOption, WindOption, DensityOption);
//		
//	int currTraj = 0;
//	double debris_state[8];	//contains state vector at time of breakup PLUS the UTC time
//	
//	int cutoff_tstep;
//	for (cutoff_tstep = 0; cutoff_tstep < num_time_steps; cutoff_tstep++) {
//		if (RocketLatLonAltStorage[currTraj][cutoff_tstep][2] > altitude_we_explode_until_km) {
//			break; } }
//	cutoff_tstep--;		//back it off one
//
////	cout << "DEBUGGING CUTOFF_TSTEP WITH HARDCODE!!!!!!\n";
////	cutoff_tstep = 9;
//	cout << "num_time_steps = " << num_time_steps <<  "     cutoff_tstep = " << cutoff_tstep << endl;
//	
//
//	int backwardsStep = 10;
//	// Go backwards so that making footprints is faster
//	for (int tstep = cutoff_tstep; tstep >= 0; tstep = tstep - backwardsStep) {
//		cout << "current debris timestep = " << tstep << endl;
//
//		// Copy the state vector into debris_state
//		memcpy((void *) debris_state, (void *) &State_Vector_Storage_Vec[currTraj][tstep][0], 7*sizeof(double));
//		debris_state[7] = Initial_UTC + RocketStateTimes[tstep]/(24*3600);	// And append the UTC time
//		
//		// Blow things up (generating a new debris profile every time)
//		for (int ix = 0; ix < num_debris_runs_per_explosion; ix++) {
//			TestCatalog.GenerateRandomPieces();
////			TestCatalog.ReduceRandomPieces(5.,0.);	//CHOPPING OUT ALL THE HIGH BALLISTIC COEFFICIENTS
//			
//			// Propagate debris
//			MyDebris.initialize_random_atmosphere();	//This isn't really the right thing to do...should have a consistent wind profile for every individual run
//															// What SHOULD be changing here are the debris properties, not the atmosphere.
//															// Though in some ways this gives us maximum variation which is nice.
//			DebrisLatLonAltStorage = MyDebris.Propagate_Debris_From_Catalog(debris_state, TestCatalog, debris_delta_t, (unsigned int) tstep);
//			
//			// Update debris variables
//			DebrisInitialUTC = debris_state[7];
//			debris_num_time_steps = MyDebris.get_num_time_steps();
//
//			//I know that I'm using a constant time step so I can generate DebrisStateTimes on my own
//			if (DebrisStateTimes == 0) {
//				delete DebrisStateTimes; }
//			DebrisStateTimes = new double [debris_num_time_steps];
//			for (int debrisDT = 0; debrisDT < debris_num_time_steps; debrisDT++) {
//				DebrisStateTimes[debrisDT] = debrisDT*debris_delta_t; }
//			DebrisFinalUTC = DebrisInitialUTC + DebrisStateTimes[debris_num_time_steps-1]/(3600.*24);
//
//
//
//			vector<vector<Point> > total_points_at = write_all_points_debris(tstepMinutes);
//
//
//			my_footprint.append_to_existing_footprint(total_points_at, DebrisInitialUTC, tstepMinutes);
//		
//		} }
//	
//
//
////	// DEBUGGING
////	double launchTimeHours = 15.;
////	double launchTimeMinutes = 45.;
////	double offsetTimeMinutes = 0.;
////	string folderName("/Users/marian/Documents/Research/GeneratedFiles/DEBUG/");
////	my_footprint.make_facet_files(folderName, launchTimeHours*60 + launchTimeMinutes, offsetTimeMinutes, tstepMinutes);
////	// END DEBUGGING
//
//	my_footprint.store_footprint_as_vector(footprintVectorFile);  //must also save the timing information found in points files
//
//	
//	
//	return;
//}

void Architecture::GoMatlab(vector<vector<Point> > &total_points_at) {
	string matlabFile("GeneratedFiles/matlabFile.m");

	vector<Point> unique_points;
	
	ofstream outFile;
	outFile.open(matlabFile.c_str(), ios::out);
	
	int timesteps_here = (INTxx) total_points_at.size();
	for (int ix = 0; ix < timesteps_here; ix++) {
		int numPointsHere = (INTxx) total_points_at[ix].size();
		
		unique_points.push_back(total_points_at[ix][0]);
		int num_unique = 1;

		for (int jx = 1; jx < numPointsHere; jx++) {
			Point currentPt = total_points_at[ix][jx];
			
			bool isDuplicate = false;
			for (int uIX = 0; uIX < num_unique; uIX++) {
				Point testPt = unique_points[uIX];
				
				double dist = sqrt(pow(currentPt.get_x() - testPt.get_x(),2) + pow(currentPt.get_y() - testPt.get_y(),2));
				if (dist < 0.1) {
					isDuplicate = true; 
					break;	} }
			
			if (!isDuplicate) {
				unique_points.push_back(currentPt); 
				num_unique++; }
		}
		
		
		outFile << "unique_points{" << ix+1 << "} = [\n";
		for (int jx = 0; jx < num_unique; jx++) {
			outFile << "  " << unique_points[jx].get_x() << "  " << unique_points[jx].get_y() << endl; } 
		outFile << "];\n\n\n";
		
		unique_points.clear();
		num_unique = 0;
	}

	outFile.close();
	
	return;
}

vector<vector<Point> > Architecture::ExplodeAtEdges(double tstepMinutes) {
// This function is to simulate the explosion of the launch vehicle at points along some envelope around the nominal.
// In an ideal world, we would propagate state vectors to arrive at malfunction turn information and use that to explode,
//  however I don't yet know what information will be available to me so I don't want to spend too much time devising a method
//  only to find out I can't actually use / validate it.
// OPERATION: Thus, this function will simply shift every piece of debris in xyz by the amount that the point on the envelope
//		is shifted from the nominal.  This is simply a rough placeholder until we can come up with something more realistic.
	

	// Need to know the timestep along nominal in which we exploded
	// WARNING: Be careful here.  Need to ensure that these UTC values are set properly before just using them
	double nominal_time = (DebrisInitialUTC - Initial_UTC)*24*3600;
	cout << "nominal time = " << std::setprecision(12) << nominal_time << endl;			// Print this out until you know you've addressed the issue
	
	double eps = 1e-5;
	int nominal_tstep = -1;
	for (int ix = 0; ix < num_time_steps; ix++) {
		if (fabs(RocketStateTimes[ix] - nominal_time) < eps) {
			nominal_tstep = ix; } }
	if (nominal_tstep == -1) {
		cout << "ERROR: Could not find nominal tstep!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
	
	
	// Every piece will have the same offset from this timestep (for each envelope point)
	int num_nominal_plus_envelope_points = (INTxx) Single_Trajectory_Tube[0].size();
	
	vector<Point> debrisOffset;
	debrisOffset.assign(num_nominal_plus_envelope_points-1, Point());
	
	// Location of nominal trajectory
	Point nominalPoint = Single_Trajectory_Tube[nominal_tstep][0];
	
	for (int ix = 1; ix < num_nominal_plus_envelope_points; ix++) {
		Point envelopePoint = Single_Trajectory_Tube[nominal_tstep][ix];
		debrisOffset[ix-1] = envelopePoint - nominalPoint;
		debrisOffset[ix-1].set_R_local(envelopePoint.get_R_local()); }
	
	
	// ---- Modeling this part off of write_all_points_debris

	// Dump those points into the all_points vector (initialized to accomodate all (desired) time steps
	//Collect all the points
	//Allocate the vector of points that we're going to write to file
	double timeInFlight = (DebrisStateTimes[debris_num_time_steps-1] - DebrisStateTimes[0])/60.;	//minutes
	int time_steps_out = (INTxx) ceil(timeInFlight/tstepMinutes);
	
	std::vector<std::vector<Point> > total_points_at;	
	total_points_at.assign(time_steps_out,std::vector<Point>());
	
	int numPieces = (INTxx) DebrisLatLonAltStorage.size();
	for (int pc = 0; pc < numPieces; pc++) {
		int numCurrTimeSteps = (INTxx) DebrisLatLonAltStorage[pc].size();
		for (int tstep = 0; tstep < numCurrTimeSteps; tstep++) {
			// Figure out the timestep index we're going to bin this into
			int binTimeStep = (INTxx) floor(DebrisStateTimes[tstep]/(60.*tstepMinutes));
			
			// I recall the lon == 0 condition fixed some problem when not all points were initialized.
			//		I don't know if this is still a problem or not.  Keep it for now.
			if (!(DebrisLatLonAltStorage[pc][tstep][1] == 0)
				&& (binTimeStep < time_steps_out)) {
				
				Point debrisPoint;
				double temp_local_R = get_local_earth_radius(DebrisLatLonAltStorage[pc][tstep][0]);
				
				debrisPoint.set_x(DebrisLatLonAltStorage[pc][tstep][1]*R_equator);
				debrisPoint.set_y(DebrisLatLonAltStorage[pc][tstep][0]*temp_local_R);
				debrisPoint.set_z(DebrisLatLonAltStorage[pc][tstep][2]);
				
				//cout << "debrisPoint: " << debrisPoint << endl;
				for (int ix = 0; ix < num_nominal_plus_envelope_points-1; ix++) {
					Point temp_pt = debrisPoint + debrisOffset[ix];
					
					temp_pt.set_R_local(temp_local_R);	//This is obviously not right, but it'll work for now.
					//cout << "      temp_pt: " << temp_pt << endl;

					total_points_at[binTimeStep].push_back(temp_pt); }

				} } }
	
	
	
	return total_points_at;
}





//~~~~~~~~~~~~~~~~~~ !!!!GRAVEYARD!!!! ~~~~~~~~~~~~~~~~~~~~~~~

// THESE TWO FUNCTIONS ARE COMMENTED OUT BECAUSE THEY DO NOT FOLLOW THE NEW OUTPUT FILE CONVENTION!!!!!!!
//void Architecture::write_all_points(string outFileName) {
//
//	//Collect all the points
//	std::vector<std::vector<Point> > total_points_at;	//trajectory and debris points all put together by timestep
//	Point temp_pt;
//
//	// Going to pretend that there's only one 10 timesteps and lump all points into alternating steps!! ~~~~~~~~~~~~~~
//	int half_num_range = 5;
//	int num_range = 2*half_num_range;
//	total_points_at.assign(2*half_num_range,std::vector<Point>());
//
//	cout << "ONLY GOING MAKE FOOTPRINTS FOR POINTS UNDER 300KM!!!!\n";
//	double cutoff = 300;
//
//	// Collecting the points from a nominal trajectory
//	int num_points_in_trajectory = RocketLatLonAltStorage[0].size();
//	for (int jx = 0; jx < half_num_range; jx++) {
//		//first load up all the nominal trajectories
//		for (int t = 0; t < num_points_in_trajectory; t++) {
//			for (int ix = 0; ix < Single_Trajectory_Tube[t].size(); ix++) {
//
//	//			temp_pt.set_x(DU*phi_over_range[run][t]);
//	//			temp_pt.set_y(DU*radius_over_range[run][t]);
//
//				if (Single_Trajectory_Tube[t][ix].get_z() < cutoff) {
//					total_points_at[2*jx].push_back(Single_Trajectory_Tube[t][ix]);			//PRETENDING HERE~~~~~~~~~~~~~~~~
//
//	//				int temp_index = total_points_at[2*jx].size();
//	//				cout << "total_pts lat = " << total_points_at[0][temp_index-1].get_y()/total_points_at[0][temp_index-1].get_R_local() * 180/PI << endl;
//				}
//
//			} } }
//
//	// Collecting the points from a single debris event
//	if (DebrisLatLonAltStorage.size() > 0) {
//	int num_debris_pieces = DebrisLatLonAltStorage[0].size();
//	cout << "numptsdebris = " << num_debris_pieces << endl;
//	for (int jx = 0; jx < half_num_range; jx++) {
//		//first load up all the nominal trajectories
//		for (int pc = 0; pc < num_debris_pieces; pc++) {
//			for (int tstep = 0; tstep < DebrisLatLonAltStorage[pc].size(); tstep++) {
////				cout << "pc = " << pc << "   size = " << DebrisLatLonAltStorage[pc].size() << endl;
//				if ((DebrisLatLonAltStorage[pc][tstep][2] < cutoff) && !(DebrisLatLonAltStorage[pc][tstep][1] == 0)) {
//					Point temp_pt;
//
//					double temp_local_R = DebrisLatLonAltStorage[pc][tstep][0];
//
//					temp_pt.set_x(DebrisLatLonAltStorage[pc][tstep][1]*R_equator);
//					temp_pt.set_y(DebrisLatLonAltStorage[pc][tstep][0]*temp_local_R);
//					temp_pt.set_z(DebrisLatLonAltStorage[pc][tstep][2]);
//					temp_pt.set_R_local(temp_local_R);
//
////					cout << "pc = " << pc << "  tstep = " << tstep << "   pt = " << temp_pt << endl;
//
//					total_points_at[2*jx].push_back(temp_pt);			//PRETENDING HERE~~~~~~~~~~~~~~~~
//				}
//
//			} } } }
//
////	//Output the first timestep for debugging puproses
////	int num_total_pts = total_points_at[0].size();
////	for (int ix = 0; ix < num_total_pts; ix++) {
////		cout << "ix = " << ix << "  pt = " << total_points_at[0][ix] << endl; }
//
//
//
//	ofstream outfile;
//	outfile.open(outFileName.c_str(),ios::out | ios::binary);
//
//	outfile.write((char *) &num_range, sizeof(num_range));
//
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here));
//	}
//
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point));
//	}
//
//	outfile.close();
//
//	return;
//}
//



//void Architecture::write_all_points_rocket(string outFileName, double tstepMinutes, double debrisCutoff_in) {
//	// DONT FORGET OUT OUTPUT THE INITIAL UTC TIME
//
//	double desiredTimeStep = tstepMinutes*60.;	//[sec]
//	double cutoff = 18.288;					//[km] Throw out points above this altitude threshold (should remove the equivalent condition in footprint files)
//	double debrisCutoff = debrisCutoff_in;	//[km] limit where we stop counting explosions in the SUA and treat them as re-entry
//
//	//Allocate the vector of points that we're going to write to file
//	std::vector<std::vector<Point> > total_points_at;	//[timestep][points]
//	Point temp_pt;
//
//	int num_range = num_time_steps;								//num_range is timesteps over the whole first stage
//	total_points_at.assign(num_range,std::vector<Point>());
//
//	// Collecting the points from a trajectory ellipse / tube
//	int num_points_in_trajectory = Single_Trajectory_Tube[0].size();		//num_points_in_trajectory is how many timesteps we're actually looking at...perhaps no difference
//	if (num_points_in_trajectory > 0) {										//checking to make sure ellipse has been allocated
//		for (int t = 0; t < num_points_in_trajectory; t++) {
//			for (int ix = 0; ix < Single_Trajectory_Tube[t].size(); ix++) {
//				if (Single_Trajectory_Tube[t][ix].get_z() < cutoff) {
//					total_points_at[t].push_back(Single_Trajectory_Tube[t][ix]);
//				} } } }
//
//	// This would be a good place to dump a batch of trajectory points into total_points_at
//	//
//	//
//	// ~~~~~~~~~~~~
//
//
//	ofstream outfile;
//	outfile.open(outFileName.c_str(),ios::out | ios::binary);
//
//	outfile.write((char *) &num_range, sizeof(num_range));
//	outfile.write((char *) &Initial_UTC, sizeof(Initial_UTC));
//
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
//
//	outfile.close();
//
//	return;
//}





///* This function appears to have been made obsolete. */
//void Architecture::stitch_together(double tstepMinutes, double debrisCutoff_in) {
//	//-- Possible inputs
//	double desiredTimeStep = tstepMinutes*60.;	//[sec]
//	double cutoff = 18.288;			//[km] Throw out points above this altitude threshold
//	//double debrisCutoff = debrisCutoff_in;	//[km] limit where we stop counting explosions in the SUA and treat them as re-entry
//    
//    
//	// These indices need to be in ascending order for the coming logic to work
//	std::vector<int> trajectories;	//Stores indices of trajectories to use
//	std::vector<int> trajectoryTubes;	//Stores indicies of trajectory tubes to use
//	std::vector<int> debrisEvents;	//Stores indices of debris generating events to use.
//	
//	//Some strings that define names of files to look for
//	
//	//-- End Possible inputs
//	
//	// Want to establish the total time range over which to collect points
//	//	This can either be specified as an input or we'll figure out the total possible range
//	
//	// First get the nominal trajectory to gather some basic info
//	read_single_trajectory_from_file("GeneratedFiles/trajectory_points.dat");
//	double stateTimeStart = Initial_UTC;
//	
//    //	for (int tstep = 0; tstep < num_time_steps; tstep++) {
//    //		if (RocketLatLonAltStorage[0][tstep][2] < debrisCutoff) {
//    //			debrisEvents.push_back(tstep); } }
//	
//	for (int tstep = 0; tstep < 41; tstep++) {
//		debrisEvents.push_back(tstep); }
//	int numDebrisFilesToRead = (INTxx) debrisEvents.size();
//	
//    
//	// Get the highest numbered debris event first to establish the (likely) longest time
//	//  Maybe even just read the debris events in reverse (pop the indices off the vector)
//	int currDebris = debrisEvents.back();
//	read_debris_from_file(currDebris);
//    
//	double debrisTimeStart = Initial_UTC;
//	double finalDebrisTimeEnd = Initial_UTC + DebrisStateTimes[debris_num_time_steps-1]/(24*3600);
//	double timeSpan = finalDebrisTimeEnd - stateTimeStart;
//	int numTimeSteps = (INTxx) ceil(timeSpan*24*3600/desiredTimeStep);
//	
//	// Resetting the debris here to make the logic cleaner for loading points into total_points
//	reset_debris();
//	
//	// Initialize relevant state times
//	double *archStateTimes = new double [numTimeSteps];
//	for (int ix = 0; ix < numTimeSteps; ix++) {
//		archStateTimes[ix] = stateTimeStart + ix*desiredTimeStep/(24*3600); }
//    
//	// Dump those points into the all_points vector (initialized to accomodate all (desired) time steps
//	//Collect all the points
//	std::vector<std::vector<Point> > total_points_at;	//trajectory and debris points all put together by timestep
//	total_points_at.assign(numTimeSteps,std::vector<Point>());
//	
//	for (int ix = 0; ix < numDebrisFilesToRead; ix++) {
//		currDebris = debrisEvents.back();
//		read_debris_from_file(currDebris);
//		cout << "debris file " << currDebris << endl;
//        
//		debrisTimeStart = Initial_UTC;
//		int numPieces = (INTxx) DebrisLatLonAltStorage.size();
//		for (int pc = 0; pc < numPieces; pc++) {
//			int numCurrTimeSteps = (INTxx) DebrisLatLonAltStorage[pc].size();
//			for (int tstep = 0; tstep < numCurrTimeSteps; tstep++) {
//				// Figure out the timestep index we're going to bin this into
//				int binTimeStep = (INTxx) floor(((debrisTimeStart-stateTimeStart)*24*3600 + DebrisStateTimes[tstep])/desiredTimeStep);
//				
//				// I recall the lon == 0 condition fixed some problem when not all points were initialized.
//				//		I don't know if this is still a problem or not.  Keep it for now.
//				if ((DebrisLatLonAltStorage[pc][tstep][2] < cutoff)
//					&& !(DebrisLatLonAltStorage[pc][tstep][1] == 0)
//					&& (binTimeStep < numTimeSteps)) {
//					
//					Point temp_pt;
//					double temp_local_R = get_local_earth_radius(DebrisLatLonAltStorage[pc][tstep][0]);
//					
//					temp_pt.set_x(DebrisLatLonAltStorage[pc][tstep][1]*R_equator);
//					temp_pt.set_y(DebrisLatLonAltStorage[pc][tstep][0]*temp_local_R);
//					temp_pt.set_z(DebrisLatLonAltStorage[pc][tstep][2]);
//					temp_pt.set_R_local(temp_local_R);
//                    
//					total_points_at[binTimeStep].push_back(temp_pt); } } }
//        
//		reset_debris();				// Those high-numbered debris files may be rather large, free up that memory
//		debrisEvents.pop_back();	//pop it so we don't read it again next time
//	}
//    
//	
//	// ---- ALL Points files should have this format -----
//	// int				number_of_time_steps
//	// double			UTC_at_start
//	// double			delta_t_in_minutes
//	// ints				for (number_of_time_steps) { number_of_points_at_this_time_step }
//	// vector<Point>	for (number_of_time_steps) { vector_of_points_at_this_time_step }
//	
//	// Write the output file used to make a footprint out of
//	ofstream outfile;
//	outfile.open("GeneratedFiles/points_to_wrap_up.dat",ios::out | ios::binary);
//	
//	outfile.write((char *) &numTimeSteps, sizeof(numTimeSteps));
//	outfile.write((char *) &stateTimeStart, sizeof(stateTimeStart));
//	outfile.write((char *) &tstepMinutes, sizeof(tstepMinutes));
//	
//	int num_points_here;
//	for (int t = 0; t < numTimeSteps; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//	
//	for (int t = 0; t < numTimeSteps; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
//	
//	outfile.close();
//	
//	
//	return;
//}



//void Architecture::make_google_earth_files(string folderName, double debrisCutoff_in){
//	// Make a plot of the nominal trajectory and the debris trajectories that represent the current stitched architecture
//
//	//double cutoff = 18.288;			//[km] Throw out points above this altitude threshold
//	double debrisCutoff = debrisCutoff_in;	//[km] limit where we stop counting explosions in the SUA and treat them as re-entry
//
//	std::vector<int> debrisEvents;	//Stores indices of debris generating events to use.
//
//	read_single_trajectory_from_file("GeneratedFiles/trajectory_points.dat");
//
//	int last_tstep = 0;
//	while ( (last_tstep < num_time_steps) & (RocketLatLonAltStorage[0][last_tstep][2] < debrisCutoff) ){
//		debrisEvents.push_back(last_tstep);
//		last_tstep++; }
//	int numDebrisFilesToRead = (INTxx) debrisEvents.size();
//
//
//	ofstream outFile;
//	outFile.open( (folderName + "GE_plot.kml").c_str() );
//
//	if (outFile.good() ) {
//		cout << "GOOGLE EARTH files = " << (folderName + "GE_plot").c_str() << endl; }
//	else {
//		cout << "GOOGLE EARTH FILE FAILED TO WRITE!!!!~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//	// Print preamble
//	string line1GE  = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
//	string line2GE  = "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
//	string line3GE  = " <Document>\n";
//	string line4GE  = "  <Style id=\"style1\">\n";
//	string line5GE  = "   <LineStyle>\n";
//	string line6GE  = "    <colorMode>random</colorMode>\n";
//	string line7GE  = "    <width>4</width>\n";
//	string line8GE  = "   </LineStyle>\n";
//	string line9GE  = "  </Style>\n";
//	string line10GE = " <name>Debris-Path</name>\n\n";
//	outFile << line1GE << line2GE << line3GE << line4GE << line5GE << line6GE << line7GE << line8GE << line9GE << line10GE;
//
//	// Print trajectory
//	char numBuf[10];
//	sprintf(numBuf, "%i", 0);
//	string numStr = numBuf;
//	string trajName = "Trajectory" + numStr;
//
//	line1GE = "  <Placemark>\n";
//	line2GE = "   <name>Piece"+trajName+"</name>\n";
//	line3GE = "   <styleUrl>#style1</styleUrl>\n";
//	line4GE = "   <MultiGeometry>\n";
//	line5GE = "    <LineString>\n";
//	line6GE = "      <altitudeMode>absolute</altitudeMode>\n";
//	line7GE = "      <coordinates>\n";
//	outFile << line1GE << line2GE << line3GE << line4GE << line5GE << line6GE << line7GE;
//
//	// GE format is Long, Lat, Alt [degrees and meters]
//	for (int tstep = 0; tstep < last_tstep; tstep++) {
//		outFile << RocketLatLonAltStorage[0][tstep][1]*180/PI << "," << RocketLatLonAltStorage[0][tstep][0]*180/PI
//			<< "," << RocketLatLonAltStorage[0][tstep][2]*1e3 << endl; }
//
//	line1GE = "    </coordinates>\n";
//	line2GE = "   </LineString>\n";
//	line3GE = "  </MultiGeometry>\n";
//	line4GE = " </Placemark>\n\n";
//	outFile << line1GE << line2GE << line3GE << line4GE;
//
//
//
//	// Put the debris loops right here
//	for (int ix = 0; ix < numDebrisFilesToRead; ix++) {
//		int currDebris = debrisEvents.back();
//		read_debris_from_file(currDebris);
//		cout << "debris file " << currDebris << endl;
//
//		int numPieces = (INTxx) DebrisLatLonAltStorage.size();
//		for (int pc = 0; pc < numPieces; pc++) {
//			int numCurrTimeSteps = (INTxx) DebrisLatLonAltStorage[pc].size();
//			// Only want to plot a fraction of these timesteps
//			int numDesiredSteps = 20;
//			int deltaStep = (INTxx) ceil(numCurrTimeSteps/(numDesiredSteps-1));
//
//
//			sprintf(numBuf, "%i_%i", currDebris,pc);
//			numStr = numBuf;
//			string debrisName = "Debris" + numStr;
//
//			line1GE = "  <Placemark>\n";
//			line2GE = "   <name>Piece"+debrisName+"</name>\n";
//			line3GE = "   <styleUrl>#style1</styleUrl>\n";
//			line4GE = "   <MultiGeometry>\n";
//			line5GE = "    <LineString>\n";
//			line6GE = "      <altitudeMode>absolute</altitudeMode>\n";
//			line7GE = "      <coordinates>\n";
//			outFile << line1GE << line2GE << line3GE << line4GE << line5GE << line6GE << line7GE;
//
//			// Output all but the final step
//			for (int tstep = 0; tstep < numDesiredSteps-1; tstep++){
//				outFile << DebrisLatLonAltStorage[pc][tstep*deltaStep][1]*180/PI << "," << DebrisLatLonAltStorage[pc][tstep*deltaStep][0]*180/PI
//				<< "," << DebrisLatLonAltStorage[pc][tstep*deltaStep][2]*1e3 << endl; }
//
//			// Make final output the very last possible tstep
//			outFile << DebrisLatLonAltStorage[pc][numCurrTimeSteps-1][1]*180/PI << "," << DebrisLatLonAltStorage[pc][numCurrTimeSteps-1][0]*180/PI
//			<< "," << DebrisLatLonAltStorage[pc][numCurrTimeSteps-1][2]*1e3 << endl;
//
//			line1GE = "    </coordinates>\n";
//			line2GE = "   </LineString>\n";
//			line3GE = "  </MultiGeometry>\n";
//			line4GE = " </Placemark>\n\n";
//			outFile << line1GE << line2GE << line3GE << line4GE;
//		}
//
//		reset_debris();				// Those high-numbered debris files may be rather large, free up that memory
//		debrisEvents.pop_back();	//pop it so we don't read it again next time
//
//	}
//
//	// End the file
//	line1GE = "  </Document>\n";
//    line2GE = "</kml>";
//	outFile << line1GE << line2GE;
//	outFile.close();
//
//	return;
//}

































//			// This section translates the state vector into something that can be easily passed off to Francisco's python code
//			if (tstep == cutoff_tstep) {
//
//				//The output should be such that you can just copy and paste into Francisco's code
//				cout << "t = " << std::setprecision(16) << cutoff_tstep << endl;
//
//				double latitude = RocketLatLonAltStorage[0][cutoff_tstep][0];
//				double longitude = RocketLatLonAltStorage[0][cutoff_tstep][1];
//				double altitude = RocketLatLonAltStorage[0][cutoff_tstep][2];
//
//				// Convert to degrees and meters
//				cout << "lat = " << latitude*180/PI << endl;
//				cout << "long = " << longitude*180/PI << endl;
//				cout << "alt = np.array(" << RocketLatLonAltStorage[0][cutoff_tstep][2]*1e3 << ")" << endl;
//
//				ColVec Veci(3,0.);
//				Veci(0) = debris_state[3];
//				Veci(1) = debris_state[4];
//				Veci(2) = debris_state[5];
//
//				//Works
//				ColVec o_R_s__ecef(3,0.);
//				o_R_s__ecef = MyDebris.Geodetic_To_ECEF(latitude*180/PI, longitude*180/PI, altitude);
////				cout << "o_R_s__ecef = " << o_R_s__ecef << endl;
//
//				MatrixM Local_C_ECEF (3,3,0.);
//				Local_C_ECEF = prod(MyDebris.Rotation(2, -latitude), MyDebris.Rotation(3, longitude));
////				cout << "Local_C_ECEF = " << Local_C_ECEF << endl;
//
//				MatrixM ECI__C__ECEF (3,3,0.);
//				ECI__C__ECEF = MyDebris.ECEF_to_ECI_rotation(DebrisInitialUTC);
////				cout << "ECI__C__ECEF = " << ECI__C__ECEF << endl;
////				cout << "debrisInitialUtc = " << DebrisInitialUTC << endl;
//
//				ColVec o_R_s__eci(3,0.);
//				o_R_s__eci = prod(ECI__C__ECEF,o_R_s__ecef);
////				cout << "o_R_s__eci = " << o_R_s__eci << endl;
//
//				ColVec eci_omega_local(3,0.);
//				eci_omega_local(2) = rotEarthRad;
////				cout << "eci_omega_local = " << eci_omega_local << endl;
//
//				MatrixM ECEF__C__Local(3,3,0.);
//				ECEF__C__Local = trans(Local_C_ECEF);
////				cout << "ECEF__C__Local = " << ECEF__C__Local << endl;
//
//
//				MatrixM ECI__C__Local;
//				ECI__C__Local = prod(ECI__C__ECEF, ECEF__C__Local);
////				cout << "ECI__C__Local = " << ECI__C__Local << endl;
//
//
//				ColVec Vlocal(3,0.);
//				Vlocal = prod(trans( ECI__C__Local),   Veci - MyDebris.Cross_Vectors(eci_omega_local, o_R_s__eci));
//
//				double Vz = Vlocal(0);
//				double Vx = Vlocal(1);
//				double Vy = Vlocal(2);
//
//				double beta = atan2(Vy, Vx);
//				double gamma = atan2(Vz, sqrt(Vx*Vx + Vy*Vy));
//
//				cout << "Vmag = " << sqrt(Vx*Vx + Vy*Vy + Vz*Vz)*1e3 << endl;
//				cout << "beta = " << beta*180/PI << endl;
//				cout << "gamma = " << gamma*180/PI << "\n\n\n";
//
//				cout << "Vx = " << Vx << endl;
//				cout << "Vy = " << Vy << endl;
//				cout << "Vz = " << Vz << endl;
//
//				cout << Veci << endl;
//				exit(2);
//			}


//			cout << "EXIT!!!" << endl;
//			exit(123);




















