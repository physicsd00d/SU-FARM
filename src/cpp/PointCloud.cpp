//
//  PointCloud.cpp
//  Prop3Dof
//
//  Created by Thomas Colvin on 9/25/13.
//  Copyright (c) 2013 Thomas Colvin. All rights reserved.
//

#include "PointCloud.h"

PointCloud::PointCloud(string pointsFile) {
	// Read in information from the pointsFile and get it ready to be used
	load_points_file_and_points(pointsFile);
    
//    // We opened the allPointsFile in load_points_file.  Now it's time to close it.
//    allPointsFile.close();
}

PointCloud::PointCloud(){
    return;
}

PointCloud::PointCloud(PointCloud *myCloud){
    
    //First load up the timing info that we can
    all_points_UTC = myCloud->all_points_UTC;
    all_points_delta_t = myCloud->all_points_delta_t;

    // Now we know how big all_points_num_range is
    all_points_num_range = myCloud->all_points_num_range;

    // Store the location info
    all_points_launchLat = myCloud->all_points_launchLat;
    all_points_launchLon = myCloud->all_points_launchLon;
    all_points_launchAzimuth = myCloud->all_points_launchAzimuth;
    
    NASkm = myCloud->NASkm;


    // Assemble the points
    
    // Ordinarily, std vector deep copies entire first layer...but here i'm not sure if i'm passing it as an object or a pointer...
    //   BE CAREFUL HERE, CHECK FOR BUG
    num_all_points_at_tstep = myCloud->num_all_points_at_tstep;  //This applies to the all_points vector, not to the footprint vector
    // Not clear that this even gets used...
    
    
//    vector<vector<Point> > all_points_total;	//[timestep][point index]
    
    // all_points_total is a nesting of std vectors and does not contain any pointers, so using an equals sign should produce a deep copy.
    all_points_total = myCloud->all_points_total;
    
    //totalNumPointsPassedInPerID is map<int, int>, so equals sign should be deep copy
    totalNumPointsPassedInPerID = myCloud->totalNumPointsPassedInPerID;

    
    
    
//    for (int tx = 0; tx < all_points_num_range; tx++){
//        int numPiecesHere = myCloud->all_points_total[tx].size();
//        for (int px = 0; px < numPiecesHere; px++) {
//            
//        }
//            
//    }
    
            
//    all_points_total = assemble_all_points_debris(flatPointArray_in, numPieces, numTimeSteps_in, maxTime, deltaTsec);

    
//    cout << "DEBUG" << endl;
//    cout << all_points_launchLat << "   " << all_points_num_range << "   " << endl;
//    cout << all_points_launchLon << "   " << all_points_UTC << "   " << endl;
//    cout << all_points_launchAzimuth << "   " << all_points_delta_t << "   " << endl;
    
//    cout << "testPt = " << all_points_total[1][10].get_gdLatDeg() << ", " << all_points_total[1][10].get_lonDeg() << endl;
    
    
    
    return;
}


PointCloud::PointCloud(vector<double> flatPointArray_in, vector<int> pointIdArray, int numPieces, vector<int> numTimeSteps_in, int maxTime, double deltaTsec,
                       double all_points_UTC_in, double all_points_delta_t_in, double timeOffsetSec,
                       double all_points_launchLat_in, double all_points_launchLon_in, double all_points_launchAzimuth_in,
                       vector<double> massArray, vector<double> areaArray, double reactionTimeSeconds, double NASkm_in){
    
//    // Cast the void arrays once here so we don't have to do it all the time later
//    double *flatPointArray = (double *) flatPointArray_in;
//    int *numTimeSteps = (int *) numTimeSteps_in;
//    int *pointIdArray = (int *) pointIdArray_in;
//    double *massArray = (double *) massArray_in;
//    double *areaArray = (double *) areaArray_in;
    
    // Expected INCOMING units
    // flatPointArray = [deg][deg][m][m/s][m/s][m/s] (lat,lon,alt,Vx,Vy,Vz)
    // massArray [kg]
    // areaArray [m^2]
    // Will eventually convert everything to km and kg
    
    NASkm = NASkm_in;
    
    //First load up the timing info that we can
    all_points_UTC = all_points_UTC_in - timeOffsetSec/(24*3600.);
    all_points_delta_t = all_points_delta_t_in; //seconds
    
    // Assemble the points
    assemble_all_points_debris(flatPointArray_in, pointIdArray, massArray, areaArray, numPieces, numTimeSteps_in, maxTime, deltaTsec, timeOffsetSec, reactionTimeSeconds);
    
    // Now we know how big all_points_num_range is
    all_points_num_range = all_points_total.size();
    
    // Store the location info
    all_points_launchLat = all_points_launchLat_in;
    all_points_launchLon = all_points_launchLon_in;
    all_points_launchAzimuth = all_points_launchAzimuth_in;
    
    //    // Debug
    //    PrintAllPoints();
    
    return;
}




PointCloud::~PointCloud(){
    
    
//    vector<int> num_all_points_at_tstep;	//This applies to the all_points vector, not to the footprint vector
//    vector<Point> all_points_single; //[point index]
//    vector<vector<Point> > all_points_total;	//[timestep][point index]
    
    Destruct_all_points();

    
    return;
}








void PointCloud::setTimingInfo(double all_points_UTC_in, double all_points_delta_t_in, int all_points_num_range_in){
    all_points_UTC = all_points_UTC_in;
    all_points_delta_t = all_points_delta_t_in;
    all_points_num_range = all_points_num_range_in;
    
    return;
    
}

void PointCloud::setLocationInfo(double all_points_launchLat_in, double all_points_launchLon_in, double all_points_launchAzimuth_in){
    all_points_launchLat = all_points_launchLat_in;
    all_points_launchLon = all_points_launchLon_in;
    all_points_launchAzimuth = all_points_launchAzimuth_in;
    
    return;
}

void PointCloud::Destruct_all_points() {
	if (all_points_single.size() > 0) {
		// all_points_single (recall this is loaded for each timestep)
		all_points_single.clear(); }
	
	// num_all_points_at_tstep
	if ( num_all_points_at_tstep.size() > 0 ) { num_all_points_at_tstep.clear(); }
	
	// close the file
	allPointsFile.close();
	
	// the points that were passed in as a vector
    
    int ptsSize = all_points_total.size();
    if (ptsSize > 0){
        for (int ix = 0; ix < ptsSize; ++ix) {
            all_points_total[ix].clear(); }
	
        all_points_total.clear();
    }
    
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

int PointCloud::load_points_file_and_points(string allPointsFileName) {
    load_points_file_not_points(allPointsFileName);
    
    // Allocate the space
    all_points_total.assign(all_points_num_range, vector<Point>());
    
    //pretty sure this should be footprint_num_range and not all_points_num_range
    for (int i = 0; i < all_points_num_range; i++) {
		// Load up the points from the current timestep into all_points
		load_points_at_timestep(i);
        
        // Assign the vector (hopefully)
        all_points_total[i] = all_points_single;
    }
    
    allPointsFile.close();
    
    
    return 0;
}



int PointCloud::load_points_file_not_points(string allPointsFileName) {
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


void PointCloud::load_points_at_timestep(int tstep) {
    // There are two ways to call this function.  Either you have:
    //   1.) Already set up a file to be read in with the function load_points_file
    //   or
    //   2.) You have loaded some points into all_points_total via the function append_to_existing_footprint
    
	// If we read in a file of points, then the allPointsFile will be open
	if (allPointsFile.is_open()) {
		if (lastTimeStep +1 == tstep) {
			all_points_single.clear();
			all_points_single.assign(num_all_points_at_tstep[tstep],Point());
			
			allPointsFile.read((char *) &all_points_single[0], all_points_single.size()*sizeof(Point));
			lastTimeStep = tstep;
		}
		else {
			cout << "ERROR!!!! YOURE ABUSING THE load_points_at_timestep FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			exit(112930);
		}
	}
	// If the file is not open, then that means we were directly passed a vector containing the points
	else {
		all_points_single.clear();
		all_points_single = all_points_total[tstep];
	}
    
	return;
}


int PointCloud::getNumPointsAtTstep(int tx){
//    return num_all_points_at_tstep[tx];
    return all_points_total[tx].size();
}

double PointCloud::getDeltaT(){
    return all_points_delta_t;
}

double PointCloud::getLaunchLat(){
    return all_points_launchLat;
}

double PointCloud::getLaunchLon(){
    return all_points_launchLon;
}

double PointCloud::getLaunchAzimuth(){
    return all_points_launchAzimuth;
}

int PointCloud::getPointsRange(){
    return all_points_total.size();
}

double PointCloud::getInitialUTC(){
    return all_points_UTC;
} 

//// This function loads the incoming points into the PointCloud and DISCARDS the pre-existing points.
////    IT DOES NOT MERGE the points together.
////    IT MERELY RECORDS some of the new data and determines if the skygrid will need to be recomputed
//// Checks that the timestep sizes are the same
//// Determines if the incoming points overstep the time bounds in either direction
//// If incoming has earlier UTC, then the current points have their UTC updated to the incoming value
////      and we record the number of t-steps difference between the two.
//// Does not appear to do anything with launch lat/lon/az information
////
//// Updates: all_points_UTC, all_points_num_range, all_points_total
////
//int PointCloud::incorporateNewPoints(vector<vector<Point> > &total_points_at, double NewInitialUTC, double NewTstepMinutes, double launchLat, double launchLon, double launchAzimuth){
//
//    // Checks that the timestep sizes are the same
//    if (NewTstepMinutes != all_points_delta_t){
//        cout << "ERROR!!  TIMESTEPS DON'T MATCH!!!  RETURNING FROM FUNCTION WITHOUT DOING ANYTHING!!!!!";
//        return -10;
//    }
//    
//    // Check the time range of the new points, may need to extend the time vector in either direction
//    int NewNumRange = total_points_at.size();
//	double all_points_final_UTC = all_points_UTC + all_points_delta_t*all_points_num_range/(60*24);
//	double NewFinalUTC = NewInitialUTC + NewTstepMinutes * NewNumRange/(60.*24.);
//    
//    // Determines if the incoming points overstep the time bounds in either direction
//    bool isLeading = (NewInitialUTC < all_points_UTC);
//	bool isTrailing = (NewFinalUTC > all_points_final_UTC);
//    
//    // As long as this is not zero, the grid will be recalculated.
//    int leadingTimeSteps = 0;
//    
//	if (isLeading || isTrailing){
//		//all_points is somehow outside the bounds of the existing footprint
//		
//		// how many leading time steps are needed as offset from existing first step
//		if (isLeading) {
//			leadingTimeSteps = (INTxx) ceil(fabs(NewInitialUTC - all_points_UTC)*24*60/all_points_delta_t);
////            all_points_UTC = all_points_UTC - leadingTimeSteps/(24.*60.);     //ceil() snaps the UTCs to one-minute intervals
//        } else {
//            leadingTimeSteps = -1;  // Negative value simply means there exist trailing timesteps (and no leading timesteps), so recalculate.
//        }
//		
////		// how many trailing time steps are needed as offset from existing final step
////		int trailingTimeSteps = 0;
////		if (isTrailing) {
////			//ceiling because no matter what, we need to add at least one timestep
//////			trailingTimeSteps = (INTxx) ceil(fabs(all_points_final_UTC - NewFinalUTC)*24*60/all_points_delta_t);    //Jan2014 I think this isn't needed anymore
//////			UTC_Final = all_points_final_UTC;   //don't actually need to save this i think
////        }
//
////        int CheckNewNumRange = leadingTimeSteps + ((INTxx) all_points_num_range) + trailingTimeSteps;
//		
//        cout << "lead = " << leadingTimeSteps << "     newnumrange = " << NewNumRange << endl;
//        cout << "(NewInitialUTC = " << NewInitialUTC << "    all_points_UTC = "  << all_points_UTC;
//        
//    }
//    
//    all_points_UTC = NewInitialUTC;
//
//    // Should update the time vector
//    // I guess we're not actually keeping track of the times in vector form
//    all_points_num_range = NewNumRange;
//
//    // Load up all the new points
//    all_points_total = total_points_at;
//    
//    // Return the offset number needed to sync old and new
//    // This is basically leadingTimeSteps
//    return leadingTimeSteps;
//}

int PointCloud::identifyYourself(){
    cout << "Base Class returns 0" << endl;
    return 0;
}

double PointCloud::getZBinHeight(){
    cout << "getZBinHeight is unset!" << endl;
    exit(5);
    
    return -1;
}

vector<vector<Point> > PointCloud::getAllPoints(){
    return all_points_total;
}

map<int,int>  PointCloud::getTotalNumPointsPassedInPerID(){
    return totalNumPointsPassedInPerID;
}

void PointCloud::PrintAllPoints(){
    
    cout << "Printing up to first five points at each timestep\n\n";
    int ptLimit = 5;
    
    int numTimeSteps = getPointsRange();
    for (int tx = 0; tx < numTimeSteps; tx++){
        int numPts = all_points_total[tx].size();
        
        // Find the min and max velocities
        double maxV = 0.;
        for (int pt = 0; pt < numPts; pt++){
            double Vx = all_points_total[tx][pt].Vx_km_s;
            double Vy = all_points_total[tx][pt].Vy_km_s;
            double Vz = all_points_total[tx][pt].Vz_km_s;
            
            double curV = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
            if (curV > maxV){
                maxV = curV;
            }
        }

        cout << "time = " << tx << endl;
        cout << "maxVel = " << maxV << endl;
        for (int pt = 0; pt < std::min(ptLimit, numPts); pt++){
            cout << all_points_total[tx][pt];
        }
        cout << "\n\n";
    }
    
//    for (it_t = all_points_total.begin(); it_t != all_points_total.end(); ++it_t){
//        int curTime = *it_t;
//        
//        cout << "time = " << curTime << endl;
//        for (it_point = all_points_total[curTime].begin(); it_point != all_points_total[curTime].end(); ++it_point){
//            if ( (it_point) == 4 ) { break; }
//            
//            cout << it_point->second; 
//        }
//        cout << "\n\n";
//    }
    
    return;
}



// Takes in an array, flatPointArray, that linearly contains all the debris trajectories in a row.  The vector is unpacked 
//  by knowing how many timesteps (numTimeSteps) are associated with each piece of debris.  If there's a reactionTime passed
//  in, then this will throw out all points that appear after the reactionTime has elapsed.  Puts all of the points into a 
//  large vector array of Point objects.  Loads everything into all_points_total which will get passed to SkyGrid for histogramming.
// I really want the SkyGrids to be self-contained in terms of normalization.  I was throwing out points
//  before and then trying to keep track of what I was throwing out, but not anymore.  All of the incoming
//  points are truncated at their last above-ground track.  Find the max time steps for each debris class
//  and expand all debris in the class to have those timesteps by padding the end with underground values.
void PointCloud::assemble_all_points_debris(vector<double> flatPointArray, vector<int> pointIdArray, vector<double> massArray,
                                            vector<double> areaArray, int numPieces, vector<int> numTimeSteps, int maxTime,
                                            double deltaTsec, double timeOffsetSec, double reactionTimeSeconds) {
    // Francisco's code outputs Lat / Lon / Alt as Deg / Deg / M
    // My Point stuff wants units of Rad / Rad / Km
    
    double tstepSeconds = all_points_delta_t;
    double m_2_km = 1e-3;

    // First count up how many debris tracks there are for each class and find the max timesteps
    for (int ix = 0; ix < pointIdArray.size(); ix++) {
        int curID = pointIdArray[ix];
        totalNumPointsPassedInPerID[curID] += 1;
        //maxStepsPerID[curID] = std::max(maxStepsPerID[curID], numTimeSteps[ix]);
    }
    // There is an issue when chaining multiple point clouds together into a single skygrid which is that
    //  for any given curID, the maximum time step seen for that ID will be different between the
    //  the different clouds.  These MUST all be identical otherwise the probabilities will get messed up.
    // HACK: Let all curIDs be propagated out to maxTime and choose maxTime to be consistent across all clouds.
    // This is a hack because 1.) it does not address the eventual issue of merging skygrids from different tfails
    //   and 2.) it seems like what should REALLY happen is that maxTime = reactionTime.
    // In the case of anticipating different tfails, maybe prepend points back to start of launch...
//    for (vector<int>::iterator idIter = pointIdArray.begin(); idIter)
    

//    // Print stuff too if you want.
//    for (std::map<int, int>::iterator it = totalNumPointsPassedInPerID.begin();
//         it !=totalNumPointsPassedInPerID.end(); ++it){
//        int curID = it->first;
//        int numThisID = it->second;
////        printf("curID %d has %d pieces\n", curID, numThisID);
////        printf("curID %d has %d pieces and max timesteps %d\n", curID, numThisID, maxStepsPerID[curID]);
//    }
    
    // Check to make sure the pointCloud timestep is not smaller than the propagated timestep
    if (tstepSeconds < deltaTsec){
        //error
        exit(20);
    }
    
    // Make these inputs!!!
    int elementsPerTimestep = 6;

    // Could avoid this hassle by using maps or lists...
    double timeInFlight = (timeOffsetSec + maxTime - 0.);	//seconds (assuming debrisStateTimes starts at zero)
	int time_steps_out = (INTxx) ceil(timeInFlight/tstepSeconds);
    
    bool useReactionTime;
    if (reactionTimeSeconds < 0) {
        useReactionTime = false; }
    else {
        useReactionTime = true; }
    
    if (useReactionTime){
        double timeInFlightCutoff = timeOffsetSec + reactionTimeSeconds;    //seconds
        if (timeInFlightCutoff < timeInFlight) {
            //update this if we're cutting off sooner
            time_steps_out = (INTxx) ceil(timeInFlightCutoff/tstepSeconds) + 1;  } } //plus 1 for converting from times to steps
    
    //    cout << "~~~~~~~~~~ Entering assemble_all_points_debris ~~~~~~~~~~~~~~~" << endl;
    //    cout << "flatPointArray = " << flatPointArray[6] << endl;
    //    cout << "numPIeces = " << numPieces << endl;
    //    cout << "numTsteps = " << numTimeSteps[1] << endl;
    //    cout << "maxTime = " << maxTime << endl;            // This is in seconds (from tfail)
    //    cout << "deltaTsec = " << deltaTsec << endl;
    //    cout << "timeOffsetSec = " << timeOffsetSec << endl;
    //
    //    cout << "tstepMinutes " << tstepMinutes << endl;
    //    cout << "time_steps_out " << time_steps_out << endl;
    //    cout << "~~~~~~~~~~ Leaving assemble_all_points_debris ~~~~~~~~~~~~~~~\n\n\n" << endl;
    
	std::vector<std::vector<Point> > total_points_at;
	total_points_at.assign(time_steps_out,std::vector<Point>());
	
	int timePointCounter = 0;   //this should wind up being equal to sum(numTimeSteps)*numPieces,  increment by 3 every time
	for (int pc = 0; pc < numPieces; pc++) {
        int curID = pointIdArray[pc];
        int numCurrTimeSteps = numTimeSteps[pc];
//        int numMaxSteps = maxStepsPerID[curID];
        int numMaxSteps = maxTime;          // TODO: Remove maxTime and just use timeStepsOut.  Should be maxTime = reaction + latency anways.
        double curMass = massArray[pc];
        double curArea = areaArray[pc] * pow(m_2_km,2);
        
		for (int tstep = 0; tstep < numMaxSteps; tstep++) {

            // Default assumption is that debris has landed and is below ground
            double gdlat = 0.;
            double lon = 0.;
            double curAlt = -1.;
            double Vx_km_s = 0.;
            double Vy_km_s = 0.;
            double Vz_km_s = 0.;
            
            // While there is data for these time steps, use it
            if (tstep < numCurrTimeSteps){
                // First parse the flatPointArray to get the information about this current data point
                gdlat = flatPointArray[timePointCounter] * PI/180.;
                lon = flatPointArray[timePointCounter + 1]  * PI/180.;
                curAlt = fmax(flatPointArray[timePointCounter+ 2] * m_2_km, 0.00001);    //Don't let last debris go negative, place 1cm above ground
                                                                                                //If it's negative here, it will be the last timestep by construction
                Vx_km_s = flatPointArray[timePointCounter + 3]  * m_2_km;
                Vy_km_s = flatPointArray[timePointCounter + 4]  * m_2_km;
                Vz_km_s = flatPointArray[timePointCounter + 5]  * m_2_km;
                
                timePointCounter += elementsPerTimestep;
            }
            
            // Find the mission-level timestep
            int binTimeStep = floor((0. + timeOffsetSec + deltaTsec*tstep)/(tstepSeconds));

            // Before throwing anything out, make a record that a point of curID was passed in.
            //      We want to know the total number of points that were the result of the monte carlo simulation.
            // totalNumPointsPassedInPerID[curID] += 1;
            
            
//            cout << "(t,lat,lon,alt_km) = " << tstep*deltaTsec << " " << gdlat*180./PI << " " << lon*180./PI << " " << curAlt << endl;
            
			// I recall the lon == 0 condition fixed some problem when not all points were initialized.
			if ((binTimeStep < time_steps_out)) {
                
                // Make sure that we're not losing the last step and that everything is as planned.
                // printf("[%d][%d] = %f\n", curID, tstep, curAlt);
                
                // I think the local radis doesn't get used anymore, now we use a fixed reference lat and lon
                double temp_local_R = 0.;
                Point temp_pt(gdlat, lon, curAlt, temp_local_R);    //TJCHERE
                
                temp_pt.set_id(curID);
                temp_pt.Vx_km_s = Vx_km_s;
                temp_pt.Vy_km_s = Vy_km_s;
                temp_pt.Vz_km_s = Vz_km_s;
                temp_pt.mass_kg = curMass;
                temp_pt.area_km2 = curArea;
                
                //printf("[%d][%d]-[%d][%d] = (%f,%f,%f)\n",pc, tstep, binTimeStep, time_steps_out, gdlat, lon, curAlt);
                
				total_points_at[binTimeStep].push_back(temp_pt);   } }
    
    }
    
    
    all_points_total = total_points_at;
    total_points_at.clear();
}


//void PointCloud::assemble_all_points_debris(void *flatPointArray_in, void *pointIdArray_in, void *massArray_in, void *areaArray_in,
//                                            int numPieces, void *numTimeSteps_in, int maxTime, double deltaTsec, double timeOffsetSec, double reactionTimeSeconds) {
//    // Francisco's code outputs Lat / Lon / Alt as Deg / Deg / M
//    // My Point stuff wants units of Rad / Rad / Km
//    
//    // double tstepMinutes timestep already saved in the PointCloud
////    double tstepMinutes = all_points_delta_t;
//    double tstepSeconds = all_points_delta_t;
//    double m_2_km = 1e-3;
// 
//    // Check to make sure the pointCloud timestep is not smaller than the propagated timestep
//    if (tstepSeconds < deltaTsec){
//        //error
//        exit(20);
//    }
//
//    // Make these inputs!!!
//    int elementsPerTimestep = 6;
//
//    // Cast the void arrays once here so we don't have to do it all the time later
//    double *flatPointArray = (double *) flatPointArray_in;
//    int *numTimeSteps = (int *) numTimeSteps_in;
//    int *pointIdArray = (int *) pointIdArray_in;
//    double *massArray = (double *) massArray_in;
//    double *areaArray = (double *) areaArray_in;
//
//    double timeInFlight = (timeOffsetSec + maxTime - 0.);	//seconds (assuming debrisStateTimes starts at zero)
//	int time_steps_out = (INTxx) ceil(timeInFlight/tstepSeconds);
//    
//    bool useReactionTime;
//    if (reactionTimeSeconds < 0) {
//        useReactionTime = false; }
//    else {
//        useReactionTime = true; }
//    
//    if (useReactionTime){
//        double timeInFlightCutoff = timeOffsetSec + reactionTimeSeconds;    //seconds
//        if (timeInFlightCutoff < timeInFlight) {
//            //update this if we're cutting off sooner
//            time_steps_out = (INTxx) ceil(timeInFlightCutoff/tstepSeconds);  } }
//    
////    cout << "~~~~~~~~~~ Entering assemble_all_points_debris ~~~~~~~~~~~~~~~" << endl;
////    cout << "flatPointArray = " << flatPointArray[6] << endl;
////    cout << "numPIeces = " << numPieces << endl;
////    cout << "numTsteps = " << numTimeSteps[1] << endl;
////    cout << "maxTime = " << maxTime << endl;            // This is in seconds (from tfail)
////    cout << "deltaTsec = " << deltaTsec << endl;
////    cout << "timeOffsetSec = " << timeOffsetSec << endl;
////
////    cout << "tstepMinutes " << tstepMinutes << endl;
////    cout << "time_steps_out " << time_steps_out << endl;
////    cout << "~~~~~~~~~~ Leaving assemble_all_points_debris ~~~~~~~~~~~~~~~\n\n\n" << endl;
//
//	std::vector<std::vector<Point> > total_points_at;
//	total_points_at.assign(time_steps_out,std::vector<Point>());
//	
//	int timePointCounter = 0;   //this should wind up being equal to sum(numTimeSteps)*numPieces,  increment by 3 every time
//	for (int pc = 0; pc < numPieces; pc++) {
//        int numCurrTimeSteps = numTimeSteps[pc];
//        int curID = pointIdArray[pc];
//        double curMass = massArray[pc];
//        double curArea = areaArray[pc] * pow(m_2_km,2);
//                
//		for (int tstep = 0; tstep < numCurrTimeSteps; tstep++) {           
//            int binTimeStep = floor((0. + timeOffsetSec + deltaTsec*tstep)/(tstepSeconds));
//			
//            double gdlat = flatPointArray[timePointCounter] * PI/180.;
//            double lon = flatPointArray[timePointCounter + 1]  * PI/180.;
//            double curAlt = fmax(flatPointArray[timePointCounter+ 2] * m_2_km, 0.00001);    //Don't let debris go negative, place 1cm above ground
//            
//            double Vx_km_s = flatPointArray[timePointCounter + 3]  * m_2_km;
//            double Vy_km_s = flatPointArray[timePointCounter + 4]  * m_2_km;
//            double Vz_km_s = flatPointArray[timePointCounter + 5]  * m_2_km;
//
//            timePointCounter += elementsPerTimestep;
//            
//            // Don't load up points that are in the air after the reaction time.  They don't need to be protected against.
//            //   If I'm doing this here, do i still need to chop later???
//            if (useReactionTime && (deltaTsec*tstep > reactionTimeSeconds)){
//                continue; }
//            
////            cout << "(t,lat,lon,alt_km) = " << tstep*deltaTsec << " " << gdlat*180./PI << " " << lon*180./PI << " " << curAlt << endl;
//            
//			// I recall the lon == 0 condition fixed some problem when not all points were initialized.
//			//		I don't know if this is still a problem or not.  Keep it for now.
//			if (!(lon == 0)
//				&& (binTimeStep < time_steps_out)) {
//                
//                // I think the local radis doesn't get used anymore, now we use a fixed reference lat and lon
//                double temp_local_R = 0.;
//                Point temp_pt(gdlat, lon, curAlt, temp_local_R);    //TJCHERE
//                
//                temp_pt.set_id(curID);                
//                temp_pt.Vx_km_s = Vx_km_s;
//                temp_pt.Vy_km_s = Vy_km_s;
//                temp_pt.Vz_km_s = Vz_km_s;
//                temp_pt.mass_kg = curMass;
//                temp_pt.area_km2 = curArea;
//                
//				total_points_at[binTimeStep].push_back(temp_pt);   } } }
//
//    
//    all_points_total = total_points_at;
//    total_points_at.clear();
//}

int PointCloud::ChopAfterSeconds(int numSeconds){
    int old_num_range = all_points_num_range;
    int new_num_range = ceil(1.*numSeconds/(all_points_delta_t));
    int numChoppedSteps = old_num_range - new_num_range;

    if (numChoppedSteps > 0){
        all_points_total.resize(new_num_range);
        all_points_num_range = new_num_range;
    }

    return all_points_num_range;
}


void PointCloud::ExportPointsToMatlab(string fileName, double deltaZkm){
    double max_z = NASkm;	//km
//    double max_z = 18.288;	//km
    double min_z = 0.0;
    double bin_size = deltaZkm;
    
    // Bin the points by ALTITUDE
    int num_bins = (INTxx) ceil((max_z - min_z)/bin_size);
    
    // Allocate the points vector
    vector<vector<vector<Point> > > pts;	//[tx][z_bin][point]

    pts.assign(all_points_num_range,vector<vector<Point> >());
	for (int t = 0; t < all_points_num_range; t++) {
		pts[t].assign(num_bins,vector<Point>()); }
    
    for (int tx = 0; tx < all_points_num_range; tx ++){
        
        // bin the points based on their z-location
        int num_current_points = (INTxx) all_points_total[tx].size();
        
        for (int i = 0; i < num_current_points; i++) {
            int bin_index = (INTxx) floor(all_points_total[tx][i].get_z()/bin_size);
            
        if ((bin_index >=0) && (bin_index < num_bins)) {
            //            cout << "current_points[" << bin_index << "].size = " << current_points[bin_index].size() << endl;
            pts[tx][bin_index].push_back(all_points_total[tx][i]); }
        else {
            // If you wound up here, that means the altitude of the current point is out of the range you specified for the bins.
        } } }
    
    
    // Set the precision and make sure output is set to scientific

    
    ofstream outfile;
	outfile.open(fileName.c_str(), ios::out);
    
    // Set the precision and make sure output is set to scientific
	outfile.precision(9);
	outfile << std::scientific;
    
    outfile << "clear all; close all; clc;\n\n";
    
    int curIX = 0;
    outfile << "vec = [\n";
    for (int tx = (all_points_num_range-1); tx < all_points_num_range; tx ++){
        int num_current_points = (INTxx) pts[tx][curIX].size();
        
        for (int ix = 0; ix < num_current_points; ix++){
            outfile << pts[tx][curIX][ix].get_x() << "  " << pts[tx][curIX][ix].get_y() << ";\n";
        } }
    outfile << "];\n\n\n\n";
    
    outfile << "scatter(vec(:,1), vec(:,2))\n";
    
    outfile.close();

//    cout.precision(9);
//	cout << std::scientific;
//    int curIX = 0;
//    cout << "quick debug...\n\n\n\n\n\nvec = [\n";
//    for (int tx = 0; tx < all_points_num_range; tx ++){
//        int num_current_points = (INTxx) pts[tx][curIX].size();
//
//        for (int ix = 0; ix < num_current_points; ix++){
//            cout << pts[tx][curIX][ix].get_x() << "  " << pts[tx][curIX][ix].get_y() << ";\n";
//        } }
//    cout << "];\n\n\n\n";
    

    

    
    return;
}




