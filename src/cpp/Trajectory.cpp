/*
 *  Trajectory.cpp
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 5/30/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

//StateDot Options
#define VACUUM 0     // wind and density option

#define SIMPLE_ATMOSPHERE 1 // wind option
#define GRAM_WINDS 2 // wind option

#define FIRST_STAGE 100
#define DEBRIS 200
#define SUBORBITAL 300
#define WHOLE_ROCKET 400
#define STAGE_DOWN 500
#define PRECOMPUTED 600

#define GRAM_DENSITY 50  //density option
#define CANTWELL_DENSITY 51		// density option
#define GRAM_DENSITY_UNCERTAIN 52		//density option


#define DEBUGGING 666

#include "Trajectory.h" 

#include "timer.h" 



Trajectory::Trajectory(){
    //This function is for reading in a pre-generated set of lat/lon/alts and outputing them to google earth
    cout << "This class now exists soley for visualizing lat/lon/alts" << endl;
    
    // Initialize these to default values (won't let you dump to google earth unless you set them to something meaningful!!!)
    launch_year = -1;
	launch_month = -1;
	launch_day = -1;
	launch_hours = -1;
	launch_minutes = -1;
    launch_seconds = -1;
    
    StateDot_Option = DEBUGGING;
    // You're not allowed to actually do anything else with this function.
    Properly_Initialized = false;
    

}

void Trajectory::loadPrecomputedFile(string inFileName, bool isDegrees, bool isAltMeters){
    cout << "Precomputed trajectory from file " << inFileName << endl;
    StateDot_Option = PRECOMPUTED;
    Get_TimeLongLatAlt(inFileName, isAltMeters, isDegrees);
}



/*!
 timeOffsetSec is often the failure time...just the time from launch.  This way all envelopes start at the same time.
 
 
 */
void Trajectory::loadDebris(vector<double> flatPointArray, vector<int> pointIdArray, vector<int> numTimeSteps, int maxTime, double deltaTsec, double timeOffsetSec, double reactionTimeMinutes){
//    cout << "Debris trajectories from file " << inFileName << endl;
    StateDot_Option = PRECOMPUTED;
    
    // Francisco's code outputs Lat / Lon / Alt as Deg / Deg / M
    // My Point stuff wants units of Rad / Rad / Km

    
//    vector<double> flatPointArray, vector<int> pointIdArray, vector<double> massArray,
//                                                vector<double> areaArray, int numPieces, vector<int> numTimeSteps, int maxTime,
//                                                double deltaTsec, double timeOffsetSec, double reactionTimeMinutes) {
    
//    cout << "\n\nCHECK UNITS!!!\n";
//    exit(-9);
    
    // Make these inputs!!!
    int elementsPerTimestep = 6;
    
    // Could avoid this hassle by using maps or lists...
    double tstepSeconds = deltaTsec;        // In case I ever want to scale the timesteps or make sure they line up with envelopes or whatever.
                                            //   This used to be in here but I took it out because it seems uneeded.  Check archives if interested.
    double timeInFlight = (timeOffsetSec + maxTime - 0.);	//seconds (assuming debrisStateTimes starts at zero)
    int time_steps_out = (INTxx) ceil(timeInFlight/tstepSeconds);
    
    bool useReactionTime;
    if (reactionTimeMinutes < 0) {
        useReactionTime = false; }
    else {
        useReactionTime = true; }
    
    if (useReactionTime){
        double timeInFlightCutoff = timeOffsetSec + reactionTimeMinutes*60.;    //seconds
        if (timeInFlightCutoff < timeInFlight) {
            //update this if we're cutting off sooner
            time_steps_out = (INTxx) ceil(timeInFlightCutoff/tstepSeconds);  } }
    
    // Allocate and initialize the state times
    PrecomputedStateTimes.assign(time_steps_out,0.);
    for (int tx = 0; tx < time_steps_out; tx++){
        PrecomputedStateTimes[tx] = tx*deltaTsec;
    }
    
    
    int numPieces = numTimeSteps.size();
    std::vector<double> zeroVec(3,-1.);             // When plotting, we'll throw out points with negative altitude.
    PrecomputedLatLonAltStorage.assign(numPieces, std::vector<std::vector<double> >());

    int timePointCounter = 0;   //this should wind up being equal to sum(numTimeSteps)*numPieces,  increment by elementsPerTimestep every time
    for (int pc = 0; pc < numPieces; pc++) {
        int numCurrTimeSteps = numTimeSteps[pc];
        int curID = pointIdArray[pc];
        
        // Allocate space for the current piece
        PrecomputedLatLonAltStorage[pc].assign(timeOffsetSec + numCurrTimeSteps, zeroVec);
        
        // Run through the tsteps FROM TIME OF FAILURE (given by timeOffsetSec)
        for (int tstep = 0; tstep < numCurrTimeSteps; tstep++) {
            
            // First parse the flatPointArray to get the information about this current data point
            double gdlat = flatPointArray[timePointCounter] * PI/180.;
            double lon = flatPointArray[timePointCounter + 1] * PI/180.;
            double curAlt = fmax(flatPointArray[timePointCounter+ 2] * 1.e-3, 0.01*1.e-3);    //Don't let debris go negative, place 1cm above ground
            
            timePointCounter += elementsPerTimestep;
            
            // Find the mission-level timestep
            int binTimeStep = floor((0. + timeOffsetSec + deltaTsec*tstep)/(tstepSeconds));
            
            if (binTimeStep < time_steps_out){
                
//                cout << time_steps_out << ",   binTimeStep = " << binTimeStep << endl;

                // This if-statement has two effects.
                //  * First, it will skip over the initial timesteps where the failure has not yet occurred (thanks to binTimeStep)
                //  * Second, it will skip over the timesteps that go beyond the reaction time, in the event that there is a reactiontime
                PrecomputedLatLonAltStorage[pc][binTimeStep][0] = gdlat;
                PrecomputedLatLonAltStorage[pc][binTimeStep][1] = lon;
                PrecomputedLatLonAltStorage[pc][binTimeStep][2] = curAlt;
            }
        }
    }
    


    
    
    
}


//Trajectory::Trajectory(string inFileName, string outFileName, bool isDegrees, bool isAltMeters){
//    //This function is for reading in a pre-generated set of lat/lon/alts and outputing them to google earth
//    cout << "This constructor exists soley for visualizing lat/lon/alts" << endl;
//    cout << "Statedot string is " << inFileName << endl;
//    
//    // The file I'm thinking of is in feet
//    //    bool isAltMeters = false;
//    //    bool isDegrees = true;
//    
//    StateDot_Option = PRECOMPUTED;
//    
//    
//    for (int ix = 0; ix < numTimeSteps; ix++){
//        PrecomputedStateTimes[ix] = *timeIT;
//        
//        // [trajectory][timestep][lat/lon/alt]
//        PrecomputedLatLonAltStorage[0][ix][0] = *latIT * angleCoeff;
//        PrecomputedLatLonAltStorage[0][ix][1] = *lonIT * angleCoeff;
//        PrecomputedLatLonAltStorage[0][ix][2] = *altIT * altCoeff/1000.;
//        
//        timeIT++;
//        latIT++;
//        lonIT++;
//        altIT++;
//    }
//    
////    Get_TimeLongLatAlt(inFileName, isAltMeters, isDegrees);
//    
//    //    char outFileBase[] = "GeneratedFiles/FaaHTHL.kml";
//    
////    write_to_google_earth_native(&outFileName[0], 1);
//    
//    // You're not allowed to actually do anything else with this function.
//    Properly_Initialized = false;
//}



void Trajectory::Get_TimeLongLatAlt(string fileName, bool isAltMeters, bool isDegrees) {
	
    double altCoeff = 0.3048;      // meters per foot
    if (isAltMeters){
        altCoeff = 1.; }           // meters per meter
    
    double angleCoeff = 1.;        // radians per radian
    if (isDegrees) {
        angleCoeff = PI/180.;      // radians per degree
    }
    
    cout << "About to import a lon/lat/alt trajectory file.  The format MUST!!! look like this:\n";
    cout << "FORMAT LATLON\n";
    cout << "# Exactly one comment line to explain things.\n";
    cout << "# Time     Lon     Lat     Alt\n";
    cout << "0          -106.    32.    101.\n";
    cout << "etc...\n\n\n";
    
	ifstream inFile;
	
	//read in the optimal trajectory file.  MUST USE c_str() or else this will fail!!!
	inFile.open(fileName.c_str(), ifstream::in);
	
	if (inFile.is_open()){
        cout << "LatLonAlt file file opened successfully!\n";
	}
	else {
		std::cerr << "file failed to open :(\n File was called: " << fileName << endl;	}
	
	char buffer[500];
	//Want to kill three lines
	for (int i=0; i < 3; i++){
		inFile.getline(buffer,500);}
    
    std::list<double> times;
    std::list<double> longitudes;
    std::list<double> latitudes;
    std::list<double> altitudes;
    
    while (!inFile.eof()){
        double tt, lon, lat, alt;
        inFile >> tt >> lon >> lat >> alt;
        inFile.ignore(256,'\n');    // ignore until new line
        
        // cout << "time = " << tt << endl;
        
        times.push_back(tt);
        longitudes.push_back(lon);
        latitudes.push_back(lat);
        altitudes.push_back(alt);
        
    }
    // cout << "Done with time" << endl;
    
    // Now we know how long it is, convert to arrays
    // Might have duplicate final point; don't care.
    int numTimeSteps = times.size();
    cout << "loaded numTimeSteps = " << numTimeSteps << endl;
    
    std::vector<double> zeroVec(3,-1.);
    PrecomputedStateTimes.assign(numTimeSteps,0.);
    PrecomputedLatLonAltStorage.assign(1, std::vector<std::vector<double> >());
    PrecomputedLatLonAltStorage[0].assign(numTimeSteps, zeroVec);
    
    std::list<double>::iterator timeIT = times.begin();
    std::list<double>::iterator lonIT = longitudes.begin();
    std::list<double>::iterator latIT = latitudes.begin();
    std::list<double>::iterator altIT = altitudes.begin();
    
    for (int ix = 0; ix < numTimeSteps; ix++){
        PrecomputedStateTimes[ix] = *timeIT;
        
        PrecomputedLatLonAltStorage[0][ix][0] = *latIT * angleCoeff;
        PrecomputedLatLonAltStorage[0][ix][1] = *lonIT * angleCoeff;
        PrecomputedLatLonAltStorage[0][ix][2] = *altIT * altCoeff/1000.;
        
        timeIT++;
        latIT++;
        lonIT++;
        altIT++;
    }
    
    // Free the lists
    times.clear();
    longitudes.clear();
    latitudes.clear();
    altitudes.clear();
    
    //    for (int ix = 0; ix < numTimeSteps; ix++){
    //        cout << "onceomre = " << PrecomputedStateTimes[ix] << endl;
    //    }
    
    
    
	
    //	//note, i'm loading these arrays backwards so that the altitude array will be monotonically increasing
    //	for (int i = NumPoints-1; i >= 0; i--) {
    //		wind_file >> WindAlt[i] >> RhoMean[i]
    //		>> WindUMean[i] >> WindVMean[i] >> WindWMean[i] >> RhoSd[i] >> WindUSd[i] >> WindVSd[i] >> WindWSd[i];
    //
    //        //		cout << i << "  wind alt = " << WindAlt[i] << endl;
    //	}
	
	inFile.close();
    return;
}


void Trajectory::setLaunchTime(int launch_year_in, int launch_month_in, int launch_day_in, int launch_hours_in, int launch_minutes_in, int launch_seconds_in){

    launch_year = launch_year_in;
	launch_month = launch_month_in;
	launch_day = launch_day_in;
	launch_hours = launch_hours_in;
	launch_minutes = launch_minutes_in;
    launch_seconds = launch_seconds_in;
    
    return;
}



int Trajectory::write_to_google_earth_native(string basename, int printThisMany) {
    return write_to_google_earth_native(basename, printThisMany, 1e10);
}


int Trajectory::write_to_google_earth_native(string basename, int printThisMany, double cutoffAltKM) {
    cout << "Valgrind says this whole function is one big memory leak" << endl;
    
    // This needs to be handled soon
    int num_time_steps = -1;
    vector <double> StateTimes;
    int j_limit = -1;       // This is the number of trajectories (debris, rocket, whatever) to be printed out
    // This loop is currently commented out
    
    std::vector<std::vector<std::vector<double> > > LatLonAltStorage;
    std::vector<std::vector<std::vector<double> > > LatLonAltStorageStageDown;
    bool isStageDown = false;
    
    
//	if (StateDot_Option == FIRST_STAGE) {
//        cout << "writing to google earth is not yet implemented for FIRST_STAGE Propagations" << endl;
//        return 10;
//		j_limit = num_per_batch;
//        //        num_time_steps = rocket_num_time_steps;
//        StateTimes = RocketStateTimes; }
//    else if (StateDot_Option == WHOLE_ROCKET) {
//        cout << "writing to google earth is not yet implemented for WHOLE_ROCKET Propagations" << endl;
//        //        return 10;
//        //		j_limit = (INTxx) State_Vector_Storage_Vec.size();
//        j_limit = printThisMany;
//        //        num_time_steps = debris_max_time_steps;
//        LatLonAltStorage = TransformToLatLonAlt(printThisMany);
//        isStageDown = true;
//        StateTimes = RocketStateTimes; }
//	else if (StateDot_Option == DEBRIS) {
//        cout << "writing to google earth is not yet implemented for DEBRIS Propagations?" << endl;
//        //        return 10;
//        // not sure what to do with printThisMany
//        LatLonAltStorage = DebrisLatLonAltStorage;
//		j_limit = num_debris_pieces;
//        num_time_steps = (INTxx) DebrisStateTimes.size();
//        StateTimes = DebrisStateTimes; }
//	else if (StateDot_Option == SUBORBITAL) {
//		j_limit = num_per_batch;
//        num_time_steps = (INTxx) State_Vector_Storage_Vec[0].size();
//        StateTimes = SuborbitalStateTimes; }
//    else if (StateDot_Option == PRECOMPUTED) {
//        j_limit = 1;
//        LatLonAltStorage = PrecomputedLatLonAltStorage;
//        isStageDown = false;
//        StateTimes = PrecomputedStateTimes;
//    }
    
    if (StateDot_Option == PRECOMPUTED) {
        j_limit = 1;
        LatLonAltStorage = PrecomputedLatLonAltStorage;
        isStageDown = false;
        StateTimes = PrecomputedStateTimes;
    } else {
        cout << "\n\n\nC++ YOU MESSED UP HARDCORE!!!\n";
        exit(-5);
    }
    
    if (launch_year < 0){
        cout << "\n\n\nC++ You must manually set the launch date before you can dump to GE!!!  EXITING.\n";
        cout << "launch_year = " << launch_year << endl << endl;
        exit(-5);
    }
    
    printThisMany = PrecomputedLatLonAltStorage.size();
    
    
    
    //Get all the kml stuff set up ----------------------------------------
    KmlFactory* factory(KmlFactory::GetFactory());
    kmldom::DocumentPtr document(factory->CreateDocument());
    document->set_description("This is TJC description right here");
    
    
    string rocketStyleString = "rocketStyleID";
    string rocketString = "rocket";
    string rocketString_hl = "rocket_hl";
    string rocketHref = "http://www.clker.com/cliparts/5/a/8/7/12375609571200265874pitr_Rocket_icon.svg.med.png";
    
    string stageDownStyleString = "stageDownStyleID";
    string stageDownString = "stageDown";
    string stageDownString_hl = "stageDown_hl";
    string stageDownHref = "http://www.clker.com/cliparts/5/a/8/7/12375609571200265874pitr_Rocket_icon.svg.med.png";
    
    kmlbase::Color32 stageDownColor;
    stageDownColor.set_color_abgr("ff3000c1");
    
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
    
    
    // StageDown Normal Pair     --------------------
    kmldom::PairPtr stageDownPair( factory->CreatePair());
    stageDownPair->set_key(0);
    stageDownPair->set_styleurl(stageDownString);
    
    // StageDown Highlighted Pair --------------------
    kmldom::PairPtr stageDownPair_hl( factory->CreatePair());
    stageDownPair_hl->set_key(1);
    stageDownPair_hl->set_styleurl(stageDownString_hl);
    
    kmldom::StyleMapPtr stageDownStyleMap(factory->CreateStyleMap());
    stageDownStyleMap->set_id(stageDownStyleString);
    stageDownStyleMap->add_pair(stageDownPair);
    stageDownStyleMap->add_pair(stageDownPair_hl);
    
    document->add_styleselector(stageDownStyleMap);
    
    
    // Define the style created in the StyleMap for ROCKET -------------------
    kmldom::IconStyleIconPtr rocketIcon( factory->CreateIconStyleIcon() );
    rocketIcon->set_href(rocketHref);
    
    kmldom::IconStylePtr rocketIconStylePtr( factory->CreateIconStyle() );
    rocketIconStylePtr->set_icon(rocketIcon);
    
    kmldom::StylePtr rocketStyle( factory->CreateStyle() );
    rocketStyle->set_id(rocketString);
    rocketStyle->set_iconstyle(rocketIconStylePtr);
    
    document->add_styleselector(rocketStyle);
    
    
    // Define the style created in the StyleMap for ROCKET_HL -------------------
    kmldom::IconStyleIconPtr rocketIcon_hl( factory->CreateIconStyleIcon() );
    rocketIcon_hl->set_href(rocketHref);
    
    kmldom::IconStylePtr rocketIconStylePtr_hl( factory->CreateIconStyle() );
    rocketIconStylePtr_hl->set_icon(rocketIcon_hl);
    
    kmldom::StylePtr rocketStyle_hl( factory->CreateStyle() );
    rocketStyle_hl->set_id(rocketString_hl);
    rocketStyle_hl->set_iconstyle(rocketIconStylePtr_hl);
    
    document->add_styleselector(rocketStyle_hl);
    
    
    
    // Define the style created in the StyleMap for STAGEDOWN -------------------
    kmldom::IconStyleIconPtr stageDownIcon( factory->CreateIconStyleIcon() );
    stageDownIcon->set_href(stageDownHref);
    
    
    kmldom::IconStylePtr stageDownIconStylePtr( factory->CreateIconStyle() );
    stageDownIconStylePtr->set_icon(stageDownIcon);
    stageDownIconStylePtr->set_color(stageDownColor);
    
    kmldom::StylePtr stageDownStyle( factory->CreateStyle() );
    stageDownStyle->set_id(stageDownString);
    stageDownStyle->set_iconstyle(stageDownIconStylePtr);
    
    document->add_styleselector(stageDownStyle);
    
    // Define the style created in the StyleMap for STAGEDOWN_HL -------------------
    kmldom::IconStyleIconPtr stageDownIcon_hl( factory->CreateIconStyleIcon() );
    stageDownIcon_hl->set_href(stageDownHref);
    
    kmldom::IconStylePtr stageDownIconStylePtr_hl( factory->CreateIconStyle() );
    stageDownIconStylePtr_hl->set_icon(stageDownIcon_hl);
    stageDownIconStylePtr_hl->set_color(stageDownColor);
    
    kmldom::StylePtr stageDownStyle_hl( factory->CreateStyle() );
    stageDownStyle_hl->set_id(stageDownString_hl);
    stageDownStyle_hl->set_iconstyle(stageDownIconStylePtr_hl);
    
    document->add_styleselector(stageDownStyle_hl);
    
    
    
    // Create the initial time structure
    time_t rawtime;
    time(&rawtime);
    
    struct tm * timeinfo;
    timeinfo = localtime ( &rawtime );  //Have to initialize the structure to the current local time then change it
    
//    int yyyy = launch_year - 1900;
//    int mm = launch_month - 1;         //May = 4 months since january
//    int dd = launch_day; //day of month
    
    timeinfo->tm_year = launch_year - 1900;     // Structure wants years SINCE 1900.
    timeinfo->tm_mon = launch_month - 1;        // Months SINCE January.  e.g. May = 4 months since January
    timeinfo->tm_mday = launch_day;
    timeinfo->tm_hour = launch_hours;
    timeinfo->tm_min = launch_minutes;
    timeinfo->tm_sec = launch_seconds;
    
    char buffer [80];
    strftime (buffer,80, "%FT%XZ", timeinfo);
    cout << "TimeOfLaunch Currently set at " << buffer << endl;
    
    //turn rawtime into the running local time
    rawtime = mktime(timeinfo);
    //    cout << "raw runnign time = " << rawtime << endl;
    
    
    char startTimeBuf[80];
    char endTimeBuf[80];
    
    // Write the geodetic coordinates of the rocket trajectories
    //	for (int j = 0; j < j_limit; j++) {
    //        int j=0;
    
    // I think i_limit should actually be the stateTimes length and then within the i loop we check to make sure we haven't violated anything.
    //        int i_limit = (INTxx) LatLonAltStorage[j].size();   //number of time steps that were propagated
    
    
    
    
    // Precompute limits of each storage vec
    double LatLonAltSize[printThisMany];
    for (int ixx = 0; ixx < printThisMany; ixx++){
        LatLonAltSize[ixx] = (INTxx) LatLonAltStorage[ixx].size(); }
    
    int i_limit = (INTxx) StateTimes.size();
    for (int i = 0; i < (i_limit-1); i++) {     // minus 1 because plotting time intervals
        
        // Get current geodetic coords
        //			ECI(0) = State_Vector_Storage_Vec[j][i][0];
        //			ECI(1) = State_Vector_Storage_Vec[j][i][1];
        //			ECI(2) = State_Vector_Storage_Vec[j][i][2];
        
        //			Geodetic = ECEF_To_Geodetic( prod(ECEF__C__ECI[i],ECI));
        
        // DO THIS IN A MINUTE!
        //Update timespan limits
        time_t startTime = rawtime + StateTimes[i];
        timeinfo = localtime ( &startTime );
        strftime (startTimeBuf,80, "%FT%XZ", timeinfo);
        
        time_t stopTime = rawtime + StateTimes[i+1];
        timeinfo = localtime ( &stopTime );
        strftime (endTimeBuf,80, "%FT%XZ", timeinfo);
        
        cout << "curStartTime " << startTimeBuf << endl;

        
        // Timespan info
        kmldom::TimeSpanPtr timespan(factory->CreateTimeSpan());
        timespan->set_begin(startTimeBuf);
        timespan->set_end(endTimeBuf);
        
        // I Think this is NOT a memory leak because, even though allocating new memory every time, I'm not losing
        //track of where it is, I'm simply storing that pointer in the folder a few lines later.
        PlacemarkPtr placemarkTemp(factory->CreatePlacemark());
        placemarkTemp->set_timeprimitive(timespan);
        
        MultiGeometryPtr multigeometryTemp(factory->CreateMultiGeometry());
        
        for (int curRun = 0; curRun < printThisMany; curRun++){
            
            if (i < LatLonAltSize[curRun]){
                // <coordinates>
                double curLat = LatLonAltStorage[curRun][i][0] * 180./PI;
                double curLon = LatLonAltStorage[curRun][i][1] * 180./PI;
//                double curAlt = (LatLonAltStorage[curRun][i][2]-local_elevation)*1.e3;
                double curAlt = (LatLonAltStorage[curRun][i][2])*1.e3;
                
                if (curAlt > 0){ 
                
                    CoordinatesPtr coordinates(factory->CreateCoordinates());
                    coordinates->add_latlngalt(curLat, curLon, curAlt);
                    
                    // <Point><coordinates>...
                    PointPtr point(factory->CreatePoint());
                    point->set_coordinates(coordinates);
                    point->set_altitudemode(kmldom::ALTITUDEMODE_RELATIVETOGROUND);
                    
                    // Put the points into a placemark
                    multigeometryTemp->add_geometry(point);
                    // printf("     %f, %f, %f\n", curLat, curLon, curAlt);
                }
            }
            
            
            //                placemarkTemp->set_geometry(point);
            
        }
        placemarkTemp->set_styleurl("#" + rocketStyleString);
        placemarkTemp->set_geometry(multigeometryTemp);
        
        
        
        // Put the placemark in the document
        document->add_feature(placemarkTemp);
    }
    
        
        
    // Create <kml> and give it <Placemark>.
    kmldom::KmlPtr kml = factory->CreateKml();
    kml->set_feature(document);
    
    // Give the xml_file a header (not sure how else to do this at the moment)
    std::string errors;
    std::string xml_output;
    kmlengine::KmlFilePtr kml_file = kmlengine::KmlFile::CreateFromImport(kml);
    kml_file->SerializeToString(&xml_output);
    
//    // Look at the output in the terminal if you want
//    cout << xml_output << endl;
    
    // Write the xml to file
    bool status = kmlbase::File::WriteStringToFile(xml_output, basename.c_str());
    if (status != 1){
        cout << "C++ GENERATING KML FILE FAILED SOMEHOW!  status = " << status << endl;
    } else {
        cout << "C++ KML File Successfully generated\n";
    }
    
	return 1;
}








//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Some general functions that I'm surprised didn't already exist

//double Trajectory::Mod(double x, double y) {
//	//Takes the modulo for a double, which cmath doesn't allow for some reason
//	
//	return (x/y - floor(x/y))*y;
//}
//
//ColVec Trajectory::Cross_Vectors(ColVec Vec1, ColVec Vec2) {
//	
//	// Set up the cross-product matrix
//    //	double temp[9] = {0,				-Vec1.at(2),		Vec1.at(1),
//    //					  Vec1.at(2),		0					-Vec1.at(0),
//    //					  -Vec1.at(1),		Vec1.at(0),			0			};
//    //
//    //	MatrixM CrossMatrix(3,3);
//    //
//    //	CrossMatrix = temp;
//    //	cout << CrossMatrix;
//	
//	ColVec Answer(3);
//	Answer(0) = Vec1(1)*Vec2(2) - Vec2(1)*Vec1(2);
//	Answer(1) = Vec1(2)*Vec2(0) - Vec2(2)*Vec1(0);
//	Answer(2) = Vec1(0)*Vec2(1) - Vec2(0)*Vec1(1);
//	
//	
//	return Answer;
//}

//double Trajectory::VecNorm(ColVec V) {
//	
//	return sqrt( V(0)*V(0) + V(1)*V(1) + V(2)*V(2));
//}
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// General Trajectory Functions
//
//
//ColVec Trajectory::Geodetic_To_ECEF(double gdlat, double lon, double alt){
//	
//    //	cout << "Earth Radius = " << DU << endl;
//	
//	// Convert from degrees to radians
//	gdlat = gdlat*PI/180;
//	lon	  = lon*PI/180;
//	
//	double Nlat = R_equator/sqrt(1-pow((ecc_Earth*sin(gdlat)),2));
//	
//	ColVec ECEF(3);// = new double[3];
//	ECEF(0) = (Nlat + alt)*cos(gdlat)*cos(lon);
//	ECEF(1) = (Nlat + alt)*cos(gdlat)*sin(lon);
//	ECEF(2) = (Nlat*(1-ecc_Earth*ecc_Earth) + alt)*sin(gdlat);
//	
//	return ECEF;
//}

//MatrixM Trajectory::ECEF_to_ECI_rotation(double UTC) {
//    //Returns transformation matrix
//	
//	double theta_g = Sidereal_Time(UTC);
//    //	double temp[9] = {cos(theta_g), -sin(theta_g), 0, sin(theta_g), cos(theta_g), 0, 0, 0, 1};
//	MatrixM ECI__C__ECEF(3,3);
//    //	cout << "Rotation formed from array\n";
//    //	ECI__C__ECEF = temp;
//    //	cout << ECI__C__ECEF;
//    //	cout << "Rotation formed from rotation function\n";
//	ECI__C__ECEF = Rotation(3,-theta_g);	//negative theta_g because ECEF must rotate BACKWARDS / clockwise to become ECI
//    //	cout << ECI__C__ECEF;
//	
//    
//	
//	return ECI__C__ECEF;
//}
//
//MatrixM Trajectory::Rotation(int axis, double angle) {
//    // BE VERY CAREFUL when using this function to make sure you're passing angle with the right sign.
//    // angle positive is a counterclockwise rotation.
//	MatrixM RotMat(3,3,0.);
//	
//	if (axis == 1) {
//		double y[9] = {	1,		0,			0,
//            0, cos(angle),	sin(angle),
//            0, -sin(angle), cos(angle)};
//		
//		RotMat(0,0) = y[0];
//		RotMat(0,1) = y[1];
//		RotMat(0,2) = y[2];
//		RotMat(1,0) = y[3];
//		RotMat(1,1) = y[4];
//		RotMat(1,2) = y[5];
//		RotMat(2,0) = y[6];
//		RotMat(2,1) = y[7];
//		RotMat(2,2) = y[8];
//	}
//	else if (axis == 2) {
//		double y[9] = {	cos(angle),	0, -sin(angle),
//            0,		1,		0,
//            sin(angle), 0, cos(angle)};
//		
//		RotMat(0,0) = y[0];
//		RotMat(0,1) = y[1];
//		RotMat(0,2) = y[2];
//		RotMat(1,0) = y[3];
//		RotMat(1,1) = y[4];
//		RotMat(1,2) = y[5];
//		RotMat(2,0) = y[6];
//		RotMat(2,1) = y[7];
//		RotMat(2,2) = y[8];
//	}
//	else if (axis == 3) {
//		double y[9] = {	cos(angle),		sin(angle), 0,
//            -sin(angle),	cos(angle), 0,
//            0,				0,		1};
//		
//		RotMat(0,0) = y[0];
//		RotMat(0,1) = y[1];
//		RotMat(0,2) = y[2];
//		RotMat(1,0) = y[3];
//		RotMat(1,1) = y[4];
//		RotMat(1,2) = y[5];
//		RotMat(2,0) = y[6];
//		RotMat(2,1) = y[7];
//		RotMat(2,2) = y[8];
//	}
//	else {
//		cout << "ERROR: You just tried to rotate about an axis that doesn't exist";
//	}
//    
//	return RotMat;
//}
//
//
//
//ColVec Trajectory::ECEF_To_Geodetic(ColVec Vec) {
//	//function [gdlat lon alt] = ECEF_To_Geodetic(x, y, z)
//    
//	//Unpack the ECI vector
//	double x = Vec(0);
//	double y = Vec(1);
//	double z = Vec(2);
//	
//	//Can pass off to other implementation from here
//	
//	double rdelta = sqrt(x*x + y*y);
//	double gdlat = asin(z/sqrt(x*x + y*y + z*z));
//    
//	// Set the tolerances
//	double tol = 1e-10;
//	double err = 50;
//	
//	double lon;
//	// I do not remember why this statement is in here.  Something to do with polar orbits it seems
//	if ((abs(x) <= tol) && (abs(y) <= tol)) {
//        //		cout << "WOAH!  If you wound up here, then you need to check that everything is okay!" << endl;
//        //		cout << Vec << endl;
//		lon = 0; }
//	else {
//		lon = atan2(y,x); }   //[deg]
//	
//    
//	// Iterate to find coordinates
//	double Nlat;
//	while( err > tol ){
//		Nlat = R_equator/sqrt(1-pow(ecc_Earth*sin(gdlat),2));
//        
//		double gdlat_new = atan((z+Nlat*sin(gdlat)*ecc_Earth*ecc_Earth)/rdelta);
//        
//		err = abs(gdlat_new - gdlat);
//		gdlat = gdlat_new;
//	}
//	
//	double alt = rdelta/cos(gdlat) - Nlat;
//	
//	//Pack up the geodetic return vector
//	ColVec Geodetic(3,1);
//	Geodetic(0) = gdlat;
//	Geodetic(1) = lon;
//	Geodetic(2) = alt;
//	
//	return Geodetic;
//}
//
//ColVec Trajectory::ECEF_To_Geodetic(double x, double y, double z) {
//	//function [gdlat lon alt] = ECEF_To_Geodetic(x, y, z)
//	
//	double rdelta = sqrt(x*x + y*y);
//	double gdlat = asin(z/sqrt(x*x + y*y + z*z));
//	
//	// Set the tolerances
//	double tol = 1e-10;
//	double err = 50;
//	
//	double lon;
//	// I do not remember why this statement is in here.  Something to do with polar orbits it seems
//	if ((abs(x) <= tol) && (abs(y) <= tol)) {
//        //		cout << "WOAH!  If you wound up here, then you need to check that everything is okay!" << endl;
//		lon = 0; }
//	else {
//		lon = atan2(y,x); }   //[deg]
//	
//	
//	// Iterate to find coordinates
//	double Nlat;
//	while( err > tol ){
//		Nlat = R_equator/sqrt(1-pow(ecc_Earth*sin(gdlat),2));
//		
//		double gdlat_new = atan((z+Nlat*sin(gdlat)*ecc_Earth*ecc_Earth)/rdelta);
//		
//		err = abs(gdlat_new - gdlat);
//		gdlat = gdlat_new;
//	}
//	
//	double alt = rdelta/cos(gdlat) - Nlat;
//	
//	//Pack up the geodetic return vector
//	ColVec Geodetic(3,1);
//	Geodetic(0) = gdlat;
//	Geodetic(1) = lon;
//	Geodetic(2) = alt;
//	
//	return Geodetic;
//}
//
//
//
//
//double Trajectory::Sidereal_Time(double UTC) {
//	//UTC is in fractional day format where midnight jan1 = 1.00
//    
//	// Greenwich Sidereal Time at 0000h 1 January 2012 UTC, from US Naval Observatory website
//	double gst2012start = 6.6706801; //% [sidereal hours] from vernal equinox
//    
//	// convert to rad
//	gst2012start = gst2012start * 2*PI/24;//    %[rad]
//    
//	// Find angle of Greenwich meridian from vernal equinox [0, 2*pi)
//	double theta_g = Mod(gst2012start + rotEarthRad*(UTC - 1)*3600*24, 2*PI);
//    
//	return theta_g;
//}
//
//double Trajectory::UTC_Time(int day, int month, int year, int hours, int minutes, int offset) {
//	
//    //	if (year != 2012){
//    //		std::cerr << "ERROR: This function is currently only good for the year 2012" << endl;
//    //	}
//	
//	int RegYear[12]  = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
//	int LeapYear[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
//	
//	int *NumDaysPerMonth;
//	if (Mod(abs(year-1600),4) == 0) {
//		NumDaysPerMonth = LeapYear;}
//	else {
//		NumDaysPerMonth = RegYear;}
//	
//	// must look at (month-1) because arrays index from zero
//	double UTC = 0;
//	for (int i = 0; i < month-1; i++) {
//		UTC = UTC + NumDaysPerMonth[i]; }
//	
//	if (day <= NumDaysPerMonth[month-1]){
//		UTC = UTC + day; }
//	else {
//		std::cerr << "ERROR: Not that many days in the month!" << endl;
//	}
//    
//	UTC = UTC + (hours+offset)/24. + minutes/(24.*60.);	
//    
//	return UTC;
//}
//
//
//
//
//
//
//
//
//
//double Trajectory::get_local_earth_radius(double gdlat /*rad*/){
//	double a = R_equator;
//	double b = R_polar;
//	
//	return sqrt((pow(a*a*cos(gdlat),2) + pow(b*b*sin(gdlat),2))/(pow(a*cos(gdlat),2) + pow(b*sin(gdlat),2)));
//}
//





//// This is just because I can't set default values with class member variables
//std::vector<std::vector<std::vector<double> > > Trajectory::TransformToLatLonAlt(int num_to_write) {
//    return TransformToLatLonAlt(num_to_write, Initial_UTC);
//}
// 
//
//
//std::vector<std::vector<std::vector<double> > > Trajectory::TransformToLatLonAlt(int num_to_write, double cur_UTC) {
//    // cur_UTC has default value of Initial_UTC
//    
//    //    int num_to_write = 1;	//would be num_per_batch if we were writing a bunch
//    
//    //    double cur_UTC = Initial_UTC;
//    std::vector< std::vector< std::vector<double> > > StateVecToTransform;	//run cases in batches of 5000
//    vector <double> StateTimes;
//    
//    if (StateDot_Option == FIRST_STAGE){
//        StateVecToTransform = State_Vector_Storage_Vec;
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == SUBORBITAL) {
//        StateVecToTransform = State_Vector_Storage_Vec;
//        StateTimes = SuborbitalStateTimes;      }
//    else if (StateDot_Option == WHOLE_ROCKET) {
//        StateVecToTransform = State_Vector_Storage_Vec;
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == STAGE_DOWN) {
//        //        cur_UTC = Initial_UTC + StageTimes[1]/(24.*3600);
//        StateVecToTransform = FallingStageStateVectorStorage;
//        StateTimes = StageDownStateTimes;   }
//    else if (StateDot_Option == DEBRIS) {
//        //        cur_UTC =
//        StateVecToTransform = Debris_State_Vector_Storage_Vec;
//        StateTimes = DebrisStateTimes;
//    }
//    
//    int StateTimesLength = (INTxx) StateTimes.size();
//    
//    // Find the transformation matrices ECEF__C__ECI at every *possible* timestep (StateTimesLength is max possible)
//	MatrixM Zero3x3 (3,3,0.);
//    std::vector< MatrixM > ECEF__C__ECI(StateTimesLength,Zero3x3);
//	for (int i = 0; i < StateTimesLength; i++) {
//		ECEF__C__ECI[i] = trans (ECEF_to_ECI_rotation(cur_UTC + StateTimes[i]/(24.*3600)) ); }
//    
//	//First generate the vector of vectors as we will want them to be read in
//	ColVec Geodetic(3,0.);
//	ColVec ECI(3,0.);
//	std::vector<double> ZeroVec3Size(3,0.);
//    std::vector<std::vector<std::vector<double> > > LatLonAltStorageTemp(num_to_write,std::vector<std::vector<double> >() );
//    
//	// Load the geodetic coordinates into the overall storage vector
//	for (int j = 0; j < num_to_write; j++) {
//		
//        //		int i_limit = (INTxx) State_Vector_Storage_Vec[j].size();
//        
//        // Only store as far as the INDIVIDUAL propagations went (should have been resized at the end of the propagation)
//		int i_limit = (INTxx) StateVecToTransform[j].size();
//		LatLonAltStorageTemp[j].assign(i_limit,ZeroVec3Size);
//		
//		for (int i = 0; i < i_limit; i++) {
//			ECI(0) = StateVecToTransform[j][i][0];
//			ECI(1) = StateVecToTransform[j][i][1];
//			ECI(2) = StateVecToTransform[j][i][2];
//			
//			Geodetic = ECEF_To_Geodetic( prod(ECEF__C__ECI[i],ECI));
//			
//			LatLonAltStorageTemp[j][i][0] = Geodetic(0);	//radians
//			LatLonAltStorageTemp[j][i][1] = Geodetic(1);	//radians
//			LatLonAltStorageTemp[j][i][2] = Geodetic(2);	//km
//			
//            //			if (i == (i_limit-1)) {
//            //				cout << "i== (i_limit-1): j = " << j << "  lat = " << Geodetic(0)*180/PI << endl; }
//            
//            if (isnan(ECI(2))){
//                cout << "hit a NAN\n\n";
//            }
//		} }
//    
//    return LatLonAltStorageTemp;
//}





/*
 *
 *
 *
 *
 ***
 *
 *
 *
 *
 */























//Trajectory::~Trajectory() {
//
////    if (Properly_Initialized) {
////        //Free the interpolators
////        switch (StateDot_Option) {
////            case FIRST_STAGE:
////                gsl_interp_accel_free (ThrustUx_acc);
////                gsl_spline_free (ThrustUx_spline);
////                gsl_interp_accel_free (ThrustUy_acc);
////                gsl_spline_free (ThrustUy_spline);
////                gsl_interp_accel_free (ThrustUz_acc);
////                gsl_spline_free (ThrustUz_spline);
////                gsl_interp_accel_free (ThrustMag_acc);
////                gsl_spline_free (ThrustMag_spline);
////                break;
////
////            case DEBRIS:
////                //Nothing yet
////                break;
////
////
////            default:
////                break;
////        }
////
////
////        gsl_interp_accel_free (RhoMean_acc);
////        gsl_spline_free (RhoMean_spline);
////        gsl_interp_accel_free (WindUMean_acc);
////        gsl_spline_free (WindUMean_spline);
////        gsl_interp_accel_free (WindVMean_acc);
////        gsl_spline_free (WindVMean_spline);
////        gsl_interp_accel_free (WindWMean_acc);
////        gsl_spline_free (WindWMean_spline);
////
////        gsl_interp_accel_free (RhoSd_acc);
////        gsl_spline_free (RhoSd_spline);
////        gsl_interp_accel_free (WindUSd_acc);
////        gsl_spline_free (WindUSd_spline);
////        gsl_interp_accel_free (WindVSd_acc);
////        gsl_spline_free (WindVSd_spline);
////        gsl_interp_accel_free (WindWSd_acc);
////        gsl_spline_free (WindWSd_spline);
////
////        gsl_interp_accel_free (Temp_acc);
////        gsl_spline_free (Temp_spline);
////
////        //done with random number generator
////        delete rand_num;
////    }
//
//	return;
//}




//double* Trajectory::TestFunction(int s) {
//	
//	double *thing = new double[s];
//	thing[0] = 12.3;
//	thing[1] = 999.001;
//	thing[2] = 54.43;
//	
//	
////	double thing[] = {12.3, 999.001, 54.43};
//	////////////
//	return thing; 
//}

//Trajectory::Trajectory(string fileName){
//    //This function is for reading in a pre-generated set of lat/lon/alts and outputing them to google earth
//    cout << "This constructor exists soley for visualizing lat/lon/alts" << endl;
//    cout << "Statedot string is " << fileName << endl;
//    
//    // The file I'm thinking of is in feet
//    bool isAltMeters = false;
//    bool isDegrees = true;
//    
//    StateDot_Option = PRECOMPUTED;
//    
//    Get_TimeLongLatAlt(fileName, isAltMeters, isDegrees);
//    
//    char outFileBase[] = "GeneratedFiles/FaaHTHL.kml";
//    write_to_google_earth_native(outFileBase, 1);
//
//    
//    
//    // You're not allowed to actually do anything else with this function.
//    Properly_Initialized = false;
//}







//Trajectory::Trajectory(char *statedot_str, char *wind_str, char *density_str) {
//
//	Properly_Initialized = false;
//	
//	// Wind and Density Options
//	char vacuum_str[] = "vacuum";
//	
//	// Wind Options
//	char simple_atm_str[] = "simple atmosphere";
//	char gram_winds_str[] = "gram winds";
//	
//	// Figure out the right wind option
//	if (!strcmp(wind_str, vacuum_str)) {
//		Wind_Option = VACUUM;
//		cout << "Wind = " << wind_str << endl; }
//	else if (!strcmp(wind_str, simple_atm_str)) {
//		Wind_Option = SIMPLE_ATMOSPHERE;
//		cout << "Wind = " << wind_str << endl; }
//	else if (!strcmp(wind_str, gram_winds_str)) {
//		Wind_Option = GRAM_WINDS;
//		cout << "Wind = " << wind_str << endl; }
//	else {
//		cout << "Invalid wind option, maybe you used capitalization?" << endl;
//	}
//
//	
//	// Density Options
//	char gram_density_str[] = "gram density";
//	char cantwell_density_str[] = "cantwell density";
//	char gram_density_uncertain_str[] = "gram density uncertain str";
//	
//	// Figure out the right density option
//	if (!strcmp(density_str, vacuum_str)) {
//		Density_Option = VACUUM;
//		cout << "Density = " << density_str << endl; }
//	
//	else if (!strcmp(density_str, gram_density_str)) {
//		Density_Option = GRAM_DENSITY;
//		cout << "Density = " << density_str << endl; }
//	
//	else if (!strcmp(density_str, cantwell_density_str)) {
//		Density_Option = CANTWELL_DENSITY;
//		cout << "Density = " << density_str << endl; }
//	
//	else if (!strcmp(density_str, gram_density_uncertain_str)) {
//		Density_Option = GRAM_DENSITY_UNCERTAIN;
//		cout << "Density = " << density_str << endl; }
//	
//	else {
//		cout << "Invalid density option, maybe you used capitalization?" << endl;
//	}
//	
//	
//	//Statedot Options
//	char first_stage_str[] = "first stage";
//    char full_rocket_str[] = "full rocket";
//	char debris_str[] = "debris";
//	char suborbital_str[] = "suborbital";
//	
//	// Figure out the right statedot option
//	if (!strcmp(statedot_str, first_stage_str)) {
//		StateDot_Option = FIRST_STAGE;
//		cout << "StateDot = " << first_stage_str << endl; }
//    else if (!strcmp(statedot_str, full_rocket_str)) {
//		StateDot_Option = WHOLE_ROCKET;
//		cout << "StateDot = " << full_rocket_str << endl; }
//	else if (!strcmp(statedot_str, debris_str)) {
//		StateDot_Option = DEBRIS;
//		cout << "StateDot = " << first_stage_str << endl; }
//	else if (!strcmp(statedot_str, suborbital_str)) {
//		StateDot_Option = SUBORBITAL;
//		cout << "StateDot = " << suborbital_str << endl; }
//	else {
//		cout << "Invalid StateDot option, maybe you used capitalization?" << endl;
//	}
//	
//	// For the moment, just going to load these in no matter what options we've choosen
//	string approx_wind_filename("Files/abridged_atmosphere.txt");
//	Get_Approximate_Wind(approx_wind_filename);
//	
//	string temperature_file("Files/Temp.txt");
//	Get_Temperature_Data(temperature_file);
//	
//	//Set up the RNG
//	rand_num = new Random_Number();
//	
//	return;
//}

//
//void Trajectory::Initialize_Stages(int num_per_batch_in, double delta_t_in, string nominal_traj_filename) {
//	
//	if (StateDot_Option != WHOLE_ROCKET) {
//		cerr << "FATAL ERROR!  Trying to initialize whole rocket, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//	
//	time_start = 0.0;
//	num_per_batch = num_per_batch_in;
//	delta_t = delta_t_in;
//	
//	Get_Nominal_Trajectory(nominal_traj_filename);
//	Aref = 5.;	//Eventually this will be read in from the nominal trajectory file.
//    Cd = 0.2;   //Eventually ...
//    cout << "WATCH OUT!!!!: Aref and Cd are hardcoded \n\n";
//	
//    // This is the number of time steps needed to capture each stage, including the zero and the terminal point
//    //  Note: Stages will overlap at the staging times
//    int StateTimesLength = 0;
//    RocketNumTimeSteps = new int [NumStages];
//    for (int i = 0; i < NumStages; i++){
//        if (i == (NumStages-1)) {
//            RocketNumTimeSteps[i] = (INTxx) ceil((ThrustTime[NumSteps-1] - StageTimes[NumStages-1])/delta_t) + 1;
//        }
//        else {
//            RocketNumTimeSteps[i] = (INTxx) ceil((StageTimes[i+1] - StageTimes[i])/delta_t) + 1; }
//        StateTimesLength += RocketNumTimeSteps[i]; }
//    StateTimesLength -= (NumStages -1);     //Correct for double-counting the interface times
//    
//    int tcount = 0;          //keeps track of the timesteps
//    int count = 0;           //keeps track of the array index
//    RocketStateTimes.assign(StateTimesLength,0.);
//    for (int i = 0; i < NumStages; i++){
////        static int tcount = 0;          //keeps track of the timesteps
////        static int count = 0;           //keeps track of the array index
//        
//        if (count >= StateTimesLength) {
//            cout << "count = " << count << "   tcount = " << tcount << endl;
//            cout << "ERROR DEBUG\n";
//        }
//
//        
//        for (int jx = 0; jx < (RocketNumTimeSteps[i]-1); jx++){
//            cout << "count = " << count << "   tcount = " << tcount << endl;
//
//            RocketStateTimes[count] = StageTimes[0] + tcount*delta_t;
//            count++;
//            tcount++; }
//        
//        
//        if (i != (NumStages-1)) {
//            RocketStateTimes[count] = StageTimes[i+1]; }
//        else {
//            // Actually we're done.  At this point, count == StateTimesLength
//            // Abandon ship.
//            //            RocketStateTimes[count] = ThrustTime[NumSteps-1];
//        }
//        count++;
//    }
//	
//	// This ONLY gets used in GRAM_DENSITY_UNCERTAIN case statements
//	//   Thus, don't (yet) need a switch statement here
//	num_uncert_density_steps = 200;
//	
//	//Allocate the storage - NOTE this is NOT a Boost vector, it's std::vector
//	std::vector<double> ZeroVecStateSize(7,0.);
//	State_Vector_Storage_Vec.assign(num_per_batch,std::vector< std::vector<double> >());
//    FallingStageStateVectorStorage.assign(num_per_batch, std::vector< std::vector<double> >());
//
////    // ---- DEBUG STUFF -----
////    FallingStageStateVectorStorage[0].assign(StateTimesLength,ZeroVecStateSize);
////
////    State_Vector_Storage_Vec[0].assign(StateTimesLength,ZeroVecStateSize);
////    // ---- DEBUG STUFF -----
//    
//	//Keep in mind that the length of each of these may get changed in the propagation (resize) later
//	for (int i = 0; i < num_per_batch; i++) {
////		State_Vector_Storage_Vec[i].assign(StateTimesLength,ZeroVecStateSize);
//		State_Vector_Storage_Vec[i].assign(StateTimesLength,std::vector<double>() );
//        for (int jx = 0; jx < StateTimesLength; jx++) {
//            State_Vector_Storage_Vec[i][jx].assign(7,0.); }
//    }
//	
//    DeltaForce.resize(3);
//
//	Properly_Initialized = true;
//	
//	return;
//}
//
//void Trajectory::Get_Random_Thrust_Offsets(int runNumber){
//    // Throw some randomness in there    
//    static int lastRun = -1;
//
//    if (runNumber == lastRun){
//        // Do nothing
//    }
//    else {
//        // Make new offset
//        DeltaForce(0) = rand() % 100;
//        DeltaForce(1) = rand() % 100;
//        DeltaForce(2) = rand() % 100;
//        double randMag = ((rand() % 100)/100.) * Fmax * (5./100.);   // Max Magnitude of disturbance is 1% of maximum thrust force (actually this looks liek 10%)
//        DeltaForce = (DeltaForce / VecNorm(DeltaForce)) * randMag;  // Turn the delta into a unit vector
//    }
//    
//    return;
//}
//
//
//
//
//
//void Trajectory::Initialize_First_Stage(int num_per_batch_in, double delta_t_in, string nominal_traj_filename) {
//	
//	if (StateDot_Option != FIRST_STAGE) {
//		cerr << "FATAL ERROR!  Trying to initialize first stage, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//	
//	time_start = 0.0;
//	num_per_batch = num_per_batch_in;
//	delta_t = delta_t_in;
//	
////	string nominal_traj_filename;
//	switch (Wind_Option) {
//		case VACUUM:
////			nominal_traj_filename = "Files/zzz_nominal_cape_vacuum.txt";
//			cout << "Running vacuum case from file = " << nominal_traj_filename << endl;
//			break;
//		case SIMPLE_ATMOSPHERE:
////			nominal_traj_filename = "Files/zzz_nominal_cape_finished.txt";
////			nominal_traj_filename = "Files/zzz_nominal_MARS.txt";
//			cout << "Running simple atmosphere case from file = " << nominal_traj_filename << endl;
//			break;
//		case GRAM_WINDS:
////			nominal_traj_filename = "Files/zzz_nominal_cape_finished.txt";
////			nominal_traj_filename = "Files/zzz_nominal_MARS.txt";
//			cout << "Running with GRAM winds atmosphere from file = " << nominal_traj_filename << endl;
//			break;	
//			
//		default:
//			cout << "ERROR: You picked an invalid statedot option!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//			break;
//	}
//	
//	Get_Nominal_Trajectory(nominal_traj_filename);
//	Aref = 5.;	//Eventually this will be read in from the nominal trajectory file.
//    Cd = 0.2;   //Eventually ...
//
//	// Determine and initialize where to propagate to
//	// Currently designed for specific case to prop up to stage separation then stop
//	int rocket_num_time_steps = (INTxx) ceil((StageTimes[1] - time_start)/delta_t) + 1;
////	RocketStateTimes = new double [rocket_num_time_steps];
//	RocketStateTimes.assign(rocket_num_time_steps,0.);
//	for (int i = 0; i < rocket_num_time_steps - 1; i++) {
//		RocketStateTimes[i] = time_start + i*delta_t; }
//	RocketStateTimes[rocket_num_time_steps - 1] = StageTimes[1];
//	
//	// This ONLY gets used in GRAM_DENSITY_UNCERTAIN case statements
//	//   Thus, don't (yet) need a switch statement here
//	num_uncert_density_steps = 200;
//	
//	//Allocate the storage - NOTE this is NOT a Boost vector, it's std::vector
//	std::vector<double> ZeroVecStateSize(7,0.);
//	State_Vector_Storage_Vec.assign(num_per_batch,std::vector< std::vector<double> >());
//	
//	for (int i = 0; i < num_per_batch; i++) {
//		State_Vector_Storage_Vec[i].assign(rocket_num_time_steps,ZeroVecStateSize); }
//	
//	Properly_Initialized = true;
//	
//	return;
//}
//
//void Trajectory::Initialize_Suborbital(double delta_t_in, string nominal_traj_filename) {
//
//	if (StateDot_Option != SUBORBITAL) {
//		cerr << "FATAL ERROR!  Trying to initialize first stage, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//
//	time_start = 0.0;
//	num_per_batch = 1;
//	delta_t = delta_t_in;
//
////	string nominal_traj_filename;
//	switch (Wind_Option) {
//		case VACUUM:
////			nominal_traj_filename = "Files/zzz_nominal_cape_vacuum.txt";
//			cout << "Running vacuum case from file = " << nominal_traj_filename << endl;
//			break;
//		case SIMPLE_ATMOSPHERE:
////			nominal_traj_filename = "Files/zzz_nominal_cape_finished.txt";
////			nominal_traj_filename = "Files/zzz_nominal_MARS.txt";
//			cout << "Running simple atmosphere case from file = " << nominal_traj_filename << endl;
//			break;
//		case GRAM_WINDS:
////			nominal_traj_filename = "Files/zzz_nominal_cape_finished.txt";
////			nominal_traj_filename = "Files/zzz_nominal_MARS.txt";
//			cout << "Running with GRAM winds atmosphere from file = " << nominal_traj_filename << endl;
//			break;
//
//		default:
//			cout << "ERROR: You picked an invalid statedot option!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//			break;
//	}
//
//	Get_Nominal_Suborbital(nominal_traj_filename);
////	Get_Nominal_Trajectory(nominal_traj_filename);
//	Aref = 5.;	//Eventually this will be read in from the nominal trajectory file.
//    Cd = 0.2;   //Eventually ...
//
////	TODO: Here
//
//	// Determine and initialize where to propagate to
//	// Currently designed for specific case to prop up to stage separation then stop
////	rocket_num_time_steps = ceil((StageTimes[1] - time_start)/delta_t) + 1;
//    
//    int rocket_num_time_steps = 1000;  //picked this to be big number, will get truncated when crash / landing occurs
////	SuborbitalStateTimes = new double [rocket_num_time_steps];
//    SuborbitalStateTimes.assign(rocket_num_time_steps, 0.);
//	for (int i = 0; i < rocket_num_time_steps - 1; i++) {
//		SuborbitalStateTimes[i] = time_start + i*delta_t; }
//	SuborbitalStateTimes[rocket_num_time_steps - 1] = StageTimes[1];
//
//	// This ONLY gets used in GRAM_DENSITY_UNCERTAIN case statements
//	//   Thus, don't (yet) need a switch statement here
//	num_uncert_density_steps = 200;
//
//	//Allocate the storage - NOTE this is NOT a Boost vector, it's std::vector
//	std::vector<double> ZeroVecStateSize(7,0.);
//	State_Vector_Storage_Vec.assign(num_per_batch,std::vector< std::vector<double> >());
//
//	for (int i = 0; i < num_per_batch; i++) {
//		State_Vector_Storage_Vec[i].assign(rocket_num_time_steps,ZeroVecStateSize); }
//
//	Properly_Initialized = true;
//
//	return;
//}
//
//
//
//
//void Trajectory::Propagate_Stage_Down_To_Ground() {
//    // Some quick error checking to make sure we're supposed to be in this function
//	if (StateDot_Option != WHOLE_ROCKET) {
//		cerr << "FATAL ERROR!  Trying to propagate some staged mass down to the ground, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//    
//    // Now we actually have to change the value to debris for the calculation to work
//    StateDot_Option = DEBRIS;
//    
//
//    
////    time_start = BreakupStateTime;       // ever using time_start not equal zero???
//    //	num_per_batch = 1;		// not currently worrying about multiple uncertain runs
//	
//	// Unsure how long it will take for debris to hit the ground
//	//  Should not take more than 24 hours!!!
//	int debris_buffer_size = (INTxx) ceil(24*3600/delta_t);
//	
//	// Vector of buffer length, where each entry is state vector at that timestep
//	std::vector<double> ZeroVecStateSize(7,0.);
//    State_Vector_Buffer_Vec.assign(debris_buffer_size, ZeroVecStateSize);
//
//    // Making the assumption that there is no uncertainty in staging time.  Change this eventually.
//    int finalIX = RocketNumTimeSteps[0]-1;
//    time_start = RocketStateTimes[finalIX];
//
//    
//    // DON'T NEED TO CALC THIS AHEAD OF TIME!  Fold into later loop
//	// Load a vector with all the possible times we're willing to consider (spans the buffer)
//	double *DebrisStateTimesBuffer = new double [debris_buffer_size];
//    StageDownStateTimes.assign(debris_buffer_size, 0.);
//	for (int i = 0; i < debris_buffer_size; i++) {
//		DebrisStateTimesBuffer[i] = time_start + i*delta_t;
//        StageDownStateTimes[i] = DebrisStateTimesBuffer[i];
//    }   //Using same delta_t as rocket
//    
//    int debris_max_time_steps = 0;
//        
//    ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
////    double temp_state[state_vec_size];
//    
//    for (int runNumber = 0; runNumber < num_per_batch; runNumber++) {
//		
////		load_density_profile(i);
////        
////		gsl_integrate_stages(state0, i);
////	}
//    
//        // Get the breakup state vectors
////        int runNumber = 0;
////        std::vector <double> BreakupStateVector = State_Vector_Storage_Vec[runNumber].back();
//        
//        // I'm imagining that the density profiles should already be loaded since this function gets called within stage integration
//        //    Initialize_density_profiles(1);		// not really worrying about this at the moment
//        //    load_density_profile(1);
//        
//        // Load up the initial state vector for propagation
////        int finalIX = (INTxx) State_Vector_Storage_Vec[runNumber].size() -1;
////        memcpy((void *) temp_state, (void *) &State_Vector_Storage_Vec[runNumber].back(), state_vec_size*sizeof(double));
//        state0(0) = State_Vector_Storage_Vec[runNumber][finalIX][0];
//        state0(1) = State_Vector_Storage_Vec[runNumber][finalIX][1];
//        state0(2) = State_Vector_Storage_Vec[runNumber][finalIX][2];
//        state0(3) = State_Vector_Storage_Vec[runNumber][finalIX][3];
//        state0(4) = State_Vector_Storage_Vec[runNumber][finalIX][4];
//        state0(5) = State_Vector_Storage_Vec[runNumber][finalIX][5];
//        state0(6) = State_Vector_Storage_Vec[runNumber][finalIX][6];
//                
//        // Returns total number of time steps that the current debris piece took to hit the ground
//        int IX = gsl_integrate_debris(state0, DebrisStateTimesBuffer);
//        
//        // Store the resulting state vectors
//        FallingStageStateVectorStorage[runNumber] = State_Vector_Buffer_Vec;
//        FallingStageStateVectorStorage[runNumber].resize(IX);
//        
////        State_Vector_Buffer_Vec.resize(IX);
////        StageDownStateTimes.resize(IX);
//        
//        // debris_max_time_steps keeps track of the maximum number of time steps used
//		if (IX > debris_max_time_steps) {
//			debris_max_time_steps = IX; }
//    }
//    
//    // Truncate state times to appropriate length
//    StageDownStateTimes.resize(debris_max_time_steps);
//    
//    cout << "(bufsize, max_time_steps) = " << debris_buffer_size << "  " << debris_max_time_steps << endl;
//    
//
//    
//    // Before returning, we must change the statedot back
//    StateDot_Option = WHOLE_ROCKET;
//        
//    //	FallingStageStateVectorStorage[run_number].assign(debris_buffer_size, ZeroVecStateSize);
//    //    FallingStageStateVectorStorage[run_number][0].assign(BreakupStateVector, BreakupStateVector + state_vec_size);
//    return;
//}
//
//
//
//std::vector<std::vector<double> > Trajectory::Propagate_Stage_Down_To_Ground(double *BreakupStateVector, double BreakupStateTime, int run_number) {
//    // Some quick error checking to make sure we're supposed to be in this function
//	if (StateDot_Option != WHOLE_ROCKET) {
//		cerr << "FATAL ERROR!  Trying to propagate some staged mass down to the ground, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//    
//    // Now we actually have to change the value to debris for the calculation to work
//    StateDot_Option = DEBRIS;
//    
//    time_start = BreakupStateTime;       // ever using time_start not equal zero???
////	num_per_batch = 1;		// not currently worrying about multiple uncertain runs
//	
//	// Unsure how long it will take for debris to hit the ground
//	//  Should not take more than 24 hours!!!
//	int debris_buffer_size = (INTxx) ceil(24*3600/delta_t);
//	
//	// Vector of buffer length, where each entry is state vector at that timestep
//	std::vector<double> ZeroVecStateSize(7,0.);
//    State_Vector_Buffer_Vec.assign(debris_buffer_size, ZeroVecStateSize);
//    
//    // DON'T NEED TO CALC THIS AHEAD OF TIME!  Fold into later loop
//	// Load a vector with all the possible times we're willing to consider (spans the buffer)
//	double *DebrisStateTimesBuffer = new double [debris_buffer_size];
//    StageDownStateTimes.assign(debris_buffer_size, 0.);
//	for (int i = 0; i < debris_buffer_size; i++) {
//		DebrisStateTimesBuffer[i] = time_start + i*delta_t;
//        StageDownStateTimes[i] = DebrisStateTimesBuffer[i];
//    }   //Using same delta_t as rocket
//    
//    // I'm imagining that the density profiles should already be loaded since this function gets called within stage integration
////    Initialize_density_profiles(1);		// not really worrying about this at the moment
////    load_density_profile(1);
//    
//    // Load up the initial state vector for propagation
//    ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
//    state0(0) = BreakupStateVector[0];
//    state0(1) = BreakupStateVector[1];
//    state0(2) = BreakupStateVector[2];
//    state0(3) = BreakupStateVector[3];
//    state0(4) = BreakupStateVector[4];
//    state0(5) = BreakupStateVector[5];
//    state0(6) = BreakupStateVector[6];
//    
//    cout << state0 << endl;
//    
//    // Returns total number of time steps that the current debris piece took to hit the ground
//    int IX = gsl_integrate_debris(state0, DebrisStateTimesBuffer);
//    State_Vector_Buffer_Vec.resize(IX);
//    StageDownStateTimes.resize(IX);
//    
////    StageStateVectorStorage = 
//    
////    std::vector<std::vector<std::vector<double> > > LatLonAltStorage(num_debris_pieces,std::vector<std::vector<double> >() );
//
//    
//    // Before returning, we must change the statedot back
//    StateDot_Option = WHOLE_ROCKET;
//
//    return State_Vector_Buffer_Vec;
//    
//    //	FallingStageStateVectorStorage[run_number].assign(debris_buffer_size, ZeroVecStateSize);
//    //    FallingStageStateVectorStorage[run_number][0].assign(BreakupStateVector, BreakupStateVector + state_vec_size);
//
//    
//}
//
//
//
//
//
//
//
///*! Initializes a single explosion instance and propagates each piece 
// 
// 
// */
//
//std::vector<std::vector<std::vector<double> > > Trajectory::Propagate_Debris_From_Catalog(double *BreakupStateVectorPlusUtcIn, Debris CurrentCatalog, double delta_t_in, unsigned int idNum) {
//// Assuming a breakup happens at BreakupStateVector, debris is generated according to CurrentCatalog.
////  Each piece of debris is propagated to the ground my marching forward in time by delta_t_in seconds.
////  idNum is most generally used to identify the specific debris event propagated, but can be used such that idNum = rocket time step index
//	
//	// Some quick error checking to make sure we're supposed to be in this function
//	if (StateDot_Option != DEBRIS) {
//		cerr << "FATAL ERROR!  Trying to initialize debris, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//	
//	// Assuming incoming state vector to be in ECI...though ECEF might be better
//	static bool needs_cleanup = false;
//	if (needs_cleanup) {
//		for (int ix = 0; ix < (int) Debris_State_Vector_Storage_Vec.size(); ix++) {
//			for (int jx = 0; jx < (int) Debris_State_Vector_Storage_Vec[ix].size(); jx++ ) {
//				Debris_State_Vector_Storage_Vec[ix].clear(); } }
//		Debris_State_Vector_Storage_Vec.clear(); 
//	}
//	needs_cleanup = true;	//so that it will get cleaned up the next time through
//	
//	// Number of pieces in current catalog
//	num_debris_pieces = CurrentCatalog.GetNumPieces();
//    
//    // SHOULD REALLY BE CALLING THE INITIALIZE DEBRIS FUNCTION FOR MOST OF THIS!!!
//    
////	cout << "Artificially reducing the size of the debris catalog!!!!!!~~~~~~~~~~~~~~~~~~~~~~" << endl;
////	num_debris_pieces = 10;
//
//	time_start = 0.0;       // ever using time_start not equal zero???
//	num_per_batch = 1;		// not currently worrying about multiple uncertain runs
//	delta_t = delta_t_in;
//	
//	// Unsure how long it will take for debris to hit the ground
//	//  Should not take more than 24 hours!!!
//	int debris_buffer_size = (INTxx) ceil(24*3600/delta_t);
//	
//	// Vector of buffer length, where each entry is state vector at that timestep
//	std::vector<double> ZeroVecStateSize(7,0.);
//	State_Vector_Buffer_Vec.assign(debris_buffer_size, ZeroVecStateSize);
//	
//	// Load a vector with all the possible times we're willing to consider (spans the buffer)
//	double *DebrisStateTimesBuffer = new double [debris_buffer_size];
//	for (int i = 0; i < debris_buffer_size; i++) {
//		DebrisStateTimesBuffer[i] = time_start + i*delta_t; }
//	
////	cout << "Time0 = " << StateTimes[0] << "   Time66 = " << StateTimes[66] << endl;
//	
//	//For debris, this will signify the number of time steps needed by the slowest of the debris pieces
//	int debris_max_time_steps = 0;
//	
//	//Allocate the storage - NOTE this is NOT a Boost vector, it's std::vector
//	Debris_State_Vector_Storage_Vec.assign(num_debris_pieces,std::vector< std::vector<double> >());  //needs cleaning
//	
//	//Create the omega vector
//	ColVec eci_omega_ecef__eci(3);
//	eci_omega_ecef__eci(0) = 0;
//	eci_omega_ecef__eci(1) = 0;	
//	eci_omega_ecef__eci(2) = rotEarthRad;	
//	
//	//Initial Conditions in ECI (All debris pieces will have same initial position vector)
//	ColVec o_R_s__eci(3,0.);	
//	ColVec eci_V_s__eci(3,0.);
//	
//	// This is just kind of being anal
//	o_R_s__eci(0) = BreakupStateVectorPlusUtcIn[0];
//	o_R_s__eci(1) = BreakupStateVectorPlusUtcIn[1];
//	o_R_s__eci(2) = BreakupStateVectorPlusUtcIn[2];
//	
//
//
//	//They will thus also have same r x omega
//	ColVec rXomega(3,0.);
//	rXomega = Cross_Vectors(o_R_s__eci, eci_omega_ecef__eci);
//	
//    cout << "BIGASS WARNING!!! ATTENUATED RANDOM OFFSETS IN DEBRIS TO SIMULATE AERODYNAMIC NON-EXPLOSIVE BREAKUP!!!!!" << endl;
//    double attenFactor = 3.;
//    
//	//~~~~~~~~~~~~~~~~ This Is the beginning of the propagation! ~~~~~~~~~~~~~~~~~~~~~~~~~~
//	DebrisPiece CurrentPiece;
//    double BreakupUTC = BreakupStateVectorPlusUtcIn[7];		//not actually launch, but breakup time in this case
//
//    
//	for (int piece_number = 0; piece_number < num_debris_pieces; piece_number++) {
//	
//		//cout << "piece " << piece_number << endl;
//		CurrentPiece = CurrentCatalog.GetDebrisPiece(piece_number);
//		
//		Aref = CurrentPiece.RefArea;
//		
//		//Because randomness is CURRENTLY isotropic, not worried too much about converting any frames
//		eci_V_s__eci(0) = BreakupStateVectorPlusUtcIn[3] + CurrentPiece.dVx/attenFactor + 0.*rXomega(0);
//		eci_V_s__eci(1) = BreakupStateVectorPlusUtcIn[4] + CurrentPiece.dVy/attenFactor + 0.*rXomega(1);
//		eci_V_s__eci(2) = BreakupStateVectorPlusUtcIn[5] + CurrentPiece.dVz/attenFactor + 0.*rXomega(2);
//		Debris_Mass0 = CurrentPiece.Mass;
//    
////		cout << "In function, mass = " << Debris_Mass0 << endl;
//		
//		// Load up the initial state vector for propagation
//		ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
//		state0(0) = o_R_s__eci(0);
//		state0(1) = o_R_s__eci(1);
//		state0(2) = o_R_s__eci(2);
//		state0(3) = eci_V_s__eci(0);
//		state0(4) = eci_V_s__eci(1);
//		state0(5) = eci_V_s__eci(2);
//		state0(6) = Debris_Mass0;
//		
//		//	cout << "initial debris vector = " << state0 << endl;
//		
//		Initialize_density_profiles(1);		// not really worrying about this at the moment
//		
//		//	for (int i = 0; i < num_per_batch; i++) {
//		
//		load_density_profile(1);
//		
//		// Returns total number of time steps that the current debris piece took to hit the ground
//		int IX = gsl_integrate_debris(state0, DebrisStateTimesBuffer);
//		
////		cout << "piece " << piece_number << " took " << IX << endl;
//		
//		// debris_max_time_steps keeps track of the maximum number of time steps used
//		if (IX > debris_max_time_steps) {
//			debris_max_time_steps = IX; }
//		
//		Debris_State_Vector_Storage_Vec[piece_number].assign(IX,ZeroVecStateSize);
//		for (int i = 0; i < IX; i++) {
//			Debris_State_Vector_Storage_Vec[piece_number][i] = State_Vector_Buffer_Vec[i]; }
//	}
//	
//	// Convert everything to LatLon
//	// Find the transformation matrices ECEF__C__ECI at every possible timestep
//    // Also construct the StateTimes vector
//	MatrixM Zero3x3 (3,3,0.);
//	std::vector< MatrixM > ECEF__C__ECI(debris_max_time_steps,Zero3x3);
//    
//    DebrisStateTimes.clear();   //Just in case
//    DebrisStateTimes.assign(debris_max_time_steps,0.);
//    
//	for (int i = 0; i < debris_max_time_steps; i++) {
//		ECEF__C__ECI[i] = trans (ECEF_to_ECI_rotation(BreakupUTC + DebrisStateTimesBuffer[i]/(24.*3600)) );
//        DebrisStateTimes[i] = DebrisStateTimesBuffer[i]; }
//	
//	
//	//First generate the vector of vectors as we will want them to be read in
//	ColVec Geodetic(3,0.);
//	ColVec ECI(3,0.);
//
//	std::vector<double> ZeroVec3Size(3,0.);
//	std::vector<std::vector<std::vector<double> > > LatLonAltStorage(num_debris_pieces,std::vector<std::vector<double> >() );
//    LatLonAltStorage = TransformToLatLonAlt(num_debris_pieces, BreakupUTC);
//
//    cout << "LatLonAlt[0][0] = " << LatLonAltStorage[0][0][0] << "   " << LatLonAltStorage[0][0][1] << "   " << LatLonAltStorage[0][0][2] << endl;
//
//    
////	// Load the geodetic coordinates into the overall storage vector
////	for (int j = 0; j < num_debris_pieces; j++) {
////		
////		int i_limit = (INTxx) Debris_State_Vector_Storage_Vec[j].size();
////		LatLonAltStorage[j].assign(i_limit,ZeroVec3Size);
////		
////		for (int i = 0; i < i_limit; i++) {
////			ECI(0) = Debris_State_Vector_Storage_Vec[j][i][0];
////			ECI(1) = Debris_State_Vector_Storage_Vec[j][i][1];
////			ECI(2) = Debris_State_Vector_Storage_Vec[j][i][2];
////			
////			Geodetic = ECEF_To_Geodetic( prod(ECEF__C__ECI[i],ECI));
////			
////			LatLonAltStorage[j][i][0] = Geodetic(0);	//radians
////			LatLonAltStorage[j][i][1] = Geodetic(1);	//radians
////			LatLonAltStorage[j][i][2] = Geodetic(2);	//km
////			
////		} }
//	
////	write_debris_to_file(idNum, LatLonAltStorage);
//
////	write_batch_to_google_earth("this field not currently used", 8);
//
//	//cout << "num tiem steps = " << num_time_steps << endl;
//    
//    //Clean up
//    delete[] DebrisStateTimesBuffer;
//	
//	
//	return LatLonAltStorage;
//}
//





//
//
//std::vector<std::vector<std::vector<double> > > Trajectory::Frisco_Propagate_Debris_From_Catalog(double *BreakupStateVectorPlusUtcIn, Debris CurrentCatalog, double delta_t_in, unsigned int idNum){
//    // Assuming a breakup happens at BreakupStateVector, debris is generated according to CurrentCatalog.
//    //  Each piece of debris is propagated to the ground my marching forward in time by delta_t_in seconds.
//    //  idNum is most generally used to identify the specific debris event propagated, but can be used such that idNum = rocket time step index
//	
//	// Some quick error checking to make sure we're supposed to be in this function
//	if (StateDot_Option != DEBRIS) {
//		cerr << "FATAL ERROR!  Trying to initialize debris, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//	
//	// Number of pieces in current catalog
//	num_debris_pieces = CurrentCatalog.GetNumPieces();
//    
//    // SHOULD REALLY BE CALLING THE INITIALIZE DEBRIS FUNCTION FOR MOST OF THIS!!!
//    
//    //	cout << "Artificially reducing the size of the debris catalog!!!!!!~~~~~~~~~~~~~~~~~~~~~~" << endl;
//    //	num_debris_pieces = 10;
//    
//	time_start = 0.0;       // ever using time_start not equal zero???
//	num_per_batch = 1;		// not currently worrying about multiple uncertain runs
//	delta_t = delta_t_in;
//	
//	// Unsure how long it will take for debris to hit the ground
//	//  Should not take more than 24 hours!!!
//	int debris_buffer_size = (INTxx) ceil(24*3600/delta_t);
//	
//	//For debris, this will signify the number of time steps needed by the slowest of the debris pieces
//	int debris_max_time_steps = 0;
//	
//	//Allocate the storage - NOTE this is NOT a Boost vector, it's std::vector
//	Debris_State_Vector_Storage_Vec.assign(num_debris_pieces,std::vector< std::vector<double> >());  //needs cleaning
//	
//	//Create the omega vector
//	ColVec eci_omega_ecef__eci(3);
//	eci_omega_ecef__eci(0) = 0;
//	eci_omega_ecef__eci(1) = 0;
//	eci_omega_ecef__eci(2) = rotEarthRad;
//	
//	//Initial Conditions in ECI (All debris pieces will have same initial position vector)
//	ColVec o_R_s__eci(3,0.);
//	ColVec eci_V_s__eci(3,0.);
//	
//	// This is just kind of being anal
//	o_R_s__eci(0) = BreakupStateVectorPlusUtcIn[0];
//	o_R_s__eci(1) = BreakupStateVectorPlusUtcIn[1];
//	o_R_s__eci(2) = BreakupStateVectorPlusUtcIn[2];
//    
//	//They will thus also have same r x omega
//	ColVec rXomega(3,0.);
//	rXomega = Cross_Vectors(o_R_s__eci, eci_omega_ecef__eci);
//	
//    cout << "BIGASS WARNING!!! ATTENUATED RANDOM OFFSETS IN DEBRIS TO SIMULATE AERODYNAMIC NON-EXPLOSIVE BREAKUP!!!!!" << endl;
//    double attenFactor = 3.;
//    
//    
//    std::vector<std::vector<std::vector<double> > > LatLonAltStorage(num_debris_pieces,std::vector<std::vector<double> >() );
//    
//	//~~~~~~~~~~~~~~~~ This Is the beginning of the propagation! ~~~~~~~~~~~~~~~~~~~~~~~~~~
//	DebrisPiece CurrentPiece;
//    double BreakupUTC = BreakupStateVectorPlusUtcIn[7];		//not actually launch, but breakup time in this case
//    
//    
//	for (int piece_number = 0; piece_number < num_debris_pieces; piece_number++) {
//        CurrentPiece = CurrentCatalog.GetDebrisPiece(piece_number);
//        cout << piece_number << " ";
//
//		Aref = CurrentPiece.RefArea;
//		
//		//Because randomness is CURRENTLY isotropic, not worried too much about converting any frames
//		eci_V_s__eci(0) = BreakupStateVectorPlusUtcIn[3] + CurrentPiece.dVx/attenFactor + 0.*rXomega(0);
//		eci_V_s__eci(1) = BreakupStateVectorPlusUtcIn[4] + CurrentPiece.dVy/attenFactor + 0.*rXomega(1);
//		eci_V_s__eci(2) = BreakupStateVectorPlusUtcIn[5] + CurrentPiece.dVz/attenFactor + 0.*rXomega(2);
//		Debris_Mass0 = CurrentPiece.Mass;
//        
//        double state0[] = { o_R_s__eci(0)*1e3, o_R_s__eci(1)*1e3, o_R_s__eci(2)*1e3, eci_V_s__eci(0)*1e3, eci_V_s__eci(1)*1e3, eci_V_s__eci(2)*1e3};
//		//        cout << "state0 xyz = " << state0[0] << "  " << state0[1] << "  " << state0[2] << endl;
//        //        cout << "state0 Vel xyz = " << state0[3] << "  " << state0[4] << "  " << state0[5] << endl;
//        //        cout << "thetag0 = " << thetag0 << endl;
//        
//        
//        
//        // Test values to be swapped out as things come online
//        double thetag0 = Sidereal_Time(BreakupUTC);
//        double debrisvel[] = {0,0,0};   //impulse velocity, already taken care of in state0
//        double mass = 100.; //Debris_Mass0;
//        double sref = 10.;  //Aref;
//        
//        int nCd = 1;
//        double *minfcd = new double[nCd];
//        minfcd[0] = 0;
//        double *cd = new double[nCd];
//        cd[0] = 1.;
//        
//        int cloption = 1;
//        
//        int nCl = 1;
//        double *minfcl = new double[nCl];
//        minfcl[0] = 1.;
//        double *cl = new double[nCl];
//        cl[0] = 0.;
//        
//        double loverd = 0.0;
//        int atmosoption = -1;
//        
//        int nList = 1;
//        
//        double *densitylist = new double[nList];
//        densitylist[0] = 0.;
//        double *ulist = new double[nList];
//        ulist[0] = 0.;
//        double *vlist = new double[nList];
//        vlist[0] = 0.;
//        double *wlist = new double[nList];
//        wlist[0] = 0.;
//        double *altitudelist = new double[nList];
//        altitudelist[0] = 0.;
//        
//        int geoptions = 0;
//        char filename[] = "debrisTrash.txt";
//        double dtinterval = delta_t;
//        int planetmodel = 1;    // Oblate Earth = 1
//        
//        
//        int ndtInterval = ceil(5*3600./dtinterval); // Since fortran will stop after 5 hours, just use that here as well
//        
//        
//        LatLonAltStorage[piece_number] = FriscoDebrisPropagation(state0,debrisvel,mass,sref,
//                 nCd,minfcd,cd,cloption,nCl,
//                 minfcl,cl,loverd,atmosoption,altitudelist,
//                 densitylist,ulist,vlist,wlist,geoptions,
//                 filename,nList, planetmodel,dtinterval,ndtInterval, thetag0);
//                
//        int IX = LatLonAltStorage[piece_number].size();
//        
////        cout << "LatLonAlt[0][0] = " << LatLonAltStorage[0][0][0] << "   " << LatLonAltStorage[0][0][1] << "   " << LatLonAltStorage[0][0][2] << endl;
//		
//        //		cout << "piece " << piece_number << " took " << IX << endl;
//		
//		// debris_max_time_steps keeps track of the maximum number of time steps used
//		if (IX > debris_max_time_steps) {
//			debris_max_time_steps = IX; }
//		
////		Debris_State_Vector_Storage_Vec[piece_number].assign(IX,ZeroVecStateSize);
////		for (int i = 0; i < IX; i++) {
////			Debris_State_Vector_Storage_Vec[piece_number][i] = State_Vector_Buffer_Vec[i]; }
//	}
//	
//	// Convert everything to LatLon
//	// Find the transformation matrices ECEF__C__ECI at every possible timestep
//    // Also construct the StateTimes vector
////	MatrixM Zero3x3 (3,3,0.);
////	std::vector< MatrixM > ECEF__C__ECI(debris_max_time_steps,Zero3x3);
//    
//    DebrisStateTimes.clear();   //Just in case
//    DebrisStateTimes.assign(debris_max_time_steps,0.);
//    
//	for (int i = 0; i < debris_max_time_steps; i++) {
////		ECEF__C__ECI[i] = trans (ECEF_to_ECI_rotation(BreakupUTC + DebrisStateTimesBuffer[i]/(24.*3600)) );
//        DebrisStateTimes[i] = time_start + i*delta_t; }
//	
//	
////	//First generate the vector of vectors as we will want them to be read in
////	ColVec Geodetic(3,0.);
////	ColVec ECI(3,0.);
////    
////	std::vector<double> ZeroVec3Size(3,0.);
////	std::vector<std::vector<std::vector<double> > > LatLonAltStorage(num_debris_pieces,std::vector<std::vector<double> >() );
////    LatLonAltStorage = TransformToLatLonAlt(num_debris_pieces, BreakupUTC);
//
//    
////    //Clean up
////    delete[] DebrisStateTimesBuffer;
//	
//	
//
//    
//    
//    
//    
////    std::vector<std::vector<std::vector<double> > > LatLonAltStorage(num_debris_pieces,std::vector<std::vector<double> >() );
//    return LatLonAltStorage;
//    
//}
//
//
//
//
//
//
//extern"C" {
//    void getspeedofsound_(double *altitude, double *speed);
//    
//    void debrispropagation_(double *finalConditions, int *numTimeSteps, double *initialState, double* DebrisVel,
//                            double *mass, double *Sref,
//                            int *nCD, double *MinfCD, double *CD,
//                            int *CLOption, int *nCL, double *MinfCL, double *CL,
//                            double *LoverD,
//                            int *atmosOption, double *altitudeList, double *densityList, double *UList, double *VList, double *WList,
//                            int *GEOptions, char *filename, int *nList, int *PlanetModel, double *dtInterval, int *ndtInterval, double *thetag0);
//    
//    
//    //    ! Input Parameters
//    //    !##############################################
//    //    double precision , dimension(3) , intent(in) :: DebrisVel
//    //    double precision , dimension(6) , intent(in) :: initialState
//    //    double precision , intent(in) :: mass,Sref, dtInterval
//    //    character*16 ,intent(in) :: filename
//    //    ! Aerodynamic inputs
//    //    integer , intent(in) :: nCd,nCL
//    //    double precision , dimension(nCD),intent(in) :: MinfCD,CD
//    //    double precision , dimension(nCL),intent(in) :: MinfCL,CL
//    //    double precision , intent(in) :: LoverD,thetag0
//    //
//    //
//    //
//    //    ! Input Options
//    //    integer , intent(in) :: GEOptions,CLOption,atmosOption,nList,PlanetModel
//    //    double precision , dimension(nList),intent(in)::altitudeList,densityList,UList,VList,WList
//    //    !##############################################
//    //
//    //    !Output Parameters
//    //    !##############################################
//    //    double precision :: latitudeFinal, longitudeFinal, altitudeFinal, VrelMag
//    //    double precision :: flag
//    //    double precision , dimension(5),intent(out) :: finalConditions
//    //
//    //    !##############################################
//    
//}
//
//std::vector<std::vector<double> > Trajectory::FriscoDebrisPropagation(double *initialState, double* DebrisVel,
//                                           double mass, double Sref,
//                                           int nCD, double *MinfCD, double *CD,
//                                           int CLOption, int nCL, double *MinfCL, double *CL,
//                                           double LoverD,
//                                           int atmosOption, double *altitudeList, double *densityList, double *UList, double *VList, double *WList,
//                                           int GEOptions, char *filename, int nList, int PlanetModel, double dtInterval, int ndtInterval, double thetag0){
//    
//    double *finalConditions;
//    finalConditions = new double [6*ndtInterval];
//    
//    int numTimeSteps = 0;
//    
//    debrispropagation_(finalConditions, &numTimeSteps, initialState,DebrisVel,&mass,&Sref,
//                       &nCD,MinfCD,CD,&CLOption,&nCL,
//                       MinfCL,CL,&LoverD,&atmosOption,altitudeList,
//                       densityList,UList,VList,WList,&GEOptions,
//                       filename,&nList, &PlanetModel,&dtInterval,&ndtInterval, &thetag0);
//    
//    // Create the return structure and load it with the points from finalConditions
//    std::vector<std::vector<double> > LatLonAlt;
//    LatLonAlt.assign(numTimeSteps, vector<double>() );
//    
//    for (int ix = 0; ix < numTimeSteps; ix++){
////        LatLonAlt[ix].assign(&finalConditions[6*ix],&finalConditions[6*ix]+3);
//        LatLonAlt[ix].assign(3,0.);
//        LatLonAlt[ix][0] = finalConditions[6*ix] * PI/180.;
//        LatLonAlt[ix][1] = finalConditions[6*ix+1] * PI/180.;
//        LatLonAlt[ix][2] = finalConditions[6*ix+2] * 1e-3;
//        
//        
////        cout << "final[" << ix << "] = " << LatLonAlt[ix][2] << endl;
//    }
//    
//    // Set the very last altitude point to zero (it's probably negative right now)
//    LatLonAlt[numTimeSteps-1][2] = 0.;
//    
//    return LatLonAlt;
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//Trajectory::Trajectory(){
//    Properly_Initialized = false;
//}
//





//
//
//void Trajectory::PropagateWholeRocket() {
//    
//    Propagate_Stages();
//    double oldCd = Cd;
//    Cd = 1.0;
//    Propagate_Stage_Down_To_Ground();
//    Cd = oldCd;
//    
//    
//    
//    return;
//}
//
//
//void Trajectory::Propagate_Stages() {
//	
//	if (Properly_Initialized == false) {
//		cerr << "FATAL ERROR: You are trying to propagate but have not initialized yet!" << endl;
//		exit(666); }
//	
//	//Rocket Initial Conditions: Stationary on the launch pad
//	ColVec a_R_s__ecef(3,0.);	//Position of rocket (s) relative to launch pad (a) is zero vector
//	ColVec ecef_V_s__ecef(3,0.);	//Velocity of rocket (s) relative to launch pad (a) is zero vector
//	
//	ColVec NoUnderscores(3,0.);	//TEST
//    
//	// Cape Canaveral coords
//	//WGS84	28° 23′ 18″ N, 80° 36′ 13″ W
////	double gdlat = 28.445455;	// geodetic degrees latitude
////	double lon = -80.564865;	// degrees longitude
////	double alt = 3e-3;			// altitude km above sea level
//	
//	ColVec o_R_a__ecef(3);
//	ColVec eci_omega_ecef__eci(3);
//    
//	o_R_a__ecef = Geodetic_To_ECEF(gdlat, lon, alt);
//	//	cout << o_R_a__ecef;
//    
//	eci_omega_ecef__eci(0) = 0;
//	eci_omega_ecef__eci(1) = 0;
//	eci_omega_ecef__eci(2) = rotEarthRad;
//	//	cout << eci_omega_ecef__eci;
//    
//	double UTC = Initial_UTC;
//	MatrixM eci__C__ecef(3,3);
//	eci__C__ecef = ECEF_to_ECI_rotation(UTC);
//	
//	
//	//Initial Conditions in ECI
//	ColVec o_R_s__eci(3,0.);
//	ColVec eci_V_s__eci(3,0.);
//	
//	o_R_s__eci = prod(eci__C__ecef,(o_R_a__ecef + a_R_s__ecef));
//	eci_V_s__eci = prod(eci__C__ecef,ecef_V_s__ecef) + Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci);
//	
//	initial_radius = VecNorm(o_R_s__eci);
//	cout << "initial radius = " << std::setprecision(16) << initial_radius << endl;
//    
//	ColVec LatLonAlt(3,0.);
//	LatLonAlt = ECEF_To_Geodetic( o_R_s__eci);
//	
//	cout << "o_R_a__ecef = \n" <<  o_R_a__ecef;
//	cout << "o_R_s__eci = \n" <<  o_R_s__eci;
//	cout << "eci_V_s__eci = \n" << eci_V_s__eci;
//	cout << "Thetag = " << Sidereal_Time(Initial_UTC) << endl;
//	cout << "LatLonAlt = " << LatLonAlt[1]*180/PI << endl;
//    
//	
//	
//	ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
//	state0(0) = o_R_s__eci(0);
//	state0(1) = o_R_s__eci(1);
//	state0(2) = o_R_s__eci(2);
//	state0(3) = eci_V_s__eci(0);
//	state0(4) = eci_V_s__eci(1);
//	state0(5) = eci_V_s__eci(2);
//	state0(6) = Mass0;
//	
//	cout << "state0 = " << state0 << endl;
//	cout << "Initial UTC = " << UTC << endl;
//    
//	// Let's just make sure that we've got an initial state vector which makes sense
//	//   Convert eci_V_s__eci into local velocity
//	cout << "lat = " << gdlat << endl;
//	cout << "lon = " << lon << endl;
//	cout << "alt = " << alt << endl;
//	cout << "initial radius = " << std::setprecision(16) << initial_radius << endl;
//    
//	// TODO: Why is this here?!  Already calculated this I think.
//	ColVec o_R_s__ecef(3,0.);
//	o_R_s__ecef = Geodetic_To_ECEF(gdlat, lon, alt);
//    //	cout << "o_R_s__ecef = " << o_R_s__ecef << endl;
//    
//	MatrixM Local_C_ECEF (3,3,0.);
//	Local_C_ECEF = prod(Rotation(2, -gdlat), Rotation(3, lon));
//    //	cout << "Local_C_ECEF = " << Local_C_ECEF << endl;
//    
//	MatrixM ECI__C__ECEF (3,3,0.);
//	ECI__C__ECEF = ECEF_to_ECI_rotation(UTC);
//    //	cout << "ECI__C__ECEF = " << ECI__C__ECEF << endl;
//    //	cout << "Initial UTC = " << UTC << endl;
//    
//    
//    
//    //	ColVec o_R_s__local(3,0.);
//    //	o_R_s__local = prod(Local_C_ECEF,o_R_s__ecef);
//    
//    //	cout << "o_R_s__eci = " << o_R_s__eci << endl;
//    
//	ColVec eci_omega_local(3,0.);
//	eci_omega_local(2) = rotEarthRad;
//    //	cout << "eci_omega_local = " << eci_omega_local << endl;
//    
//    
//    
//	MatrixM ECEF__C__Local(3,3,0.);
//	ECEF__C__Local = trans(Local_C_ECEF);
//    //	cout << "ECEF__C__Local = " << ECEF__C__Local << endl;
//    
//	MatrixM ECI__C__Local;
//	ECI__C__Local = prod(ECI__C__ECEF, ECEF__C__Local);
//    //	cout << "ECI__C__Local = " << ECI__C__Local << endl;
//    
//	ColVec Vlocal(3,1.);
//	Vlocal = prod(trans( ECI__C__Local),   eci_V_s__eci - Cross_Vectors(eci_omega_local, o_R_s__eci));
//    //	cout << "Vlocal = " << Vlocal << endl;
//    
//    
//    //	ColVec Vecef(3,0.);
//    //	Vecef = prod(ECEF__C__ECI,eci_V_s__eci);
//    //
//    //
//    //
//    //	ColVec Vlocal(3,0.);
//    //	Vlocal = prod(Local_C_ECEF,Vecef);
//    
//	//-------------- End Debugging Stuff
//    
//    
//	
//	Initialize_density_profiles(num_per_batch);
//    
//    //	timer stopwatch;
//	
//    //	stopwatch.start();
//    
//	for (int i = 0; i < num_per_batch; i++) {
//		
//		load_density_profile(i);
//        
//        Get_Random_Thrust_Offsets(i);
//        
//		gsl_integrate_stages(state0, i);
//	}
//    
//
//    // Calculate and store the azimuth of the launch!
//    // Convert the nominal trajectory (#1) to LatLon
////    std::vector<std::vector<std::vector<double> > > RocketLatLonAltStorage = TransformToLatLonAlt(1);
//    std::vector<std::vector<std::vector<double> > > RocketLatLonAltStorage = TransformToLatLonAlt(num_per_batch);
//    
//    double timeOfStaging = ceil((StageTimes[1] - StageTimes[0])/delta_t);
//    cout << "timeOfStaging = " << timeOfStaging << endl;
//	//Want semi-minor axis along direction of lateral motion, so find direction of motion
//	double LatInitial = RocketLatLonAltStorage[0][0][0];
//	double LonInitial = RocketLatLonAltStorage[0][0][1];
//	double LatFinal = RocketLatLonAltStorage[0][timeOfStaging][0];    //Using the location at stage separation
//	double LonFinal = RocketLatLonAltStorage[0][timeOfStaging][1];
//
//    Point ptInitial(LatInitial, LonInitial, 0, 0);
//    Point ptFinal(LatFinal, LonFinal, 0, 0);
//    
//    launchAzimuth = ptInitial.calcAzimuth(ptFinal);
//    cout << "azimuth = " << launchAzimuth << endl;
//    
//    // Find out when's the last time the rocket leaves the cutoff
//    int longestTime = 0.;   //timesteps, not seconds
//    double cutoff = 100;    //height of cutoff in km
//    
////    for (int px = 0; px < numPieces; px++) {
////        int numTSteps = LatLonAltStorage[px].size();
////        
////        for (int ix = 0; ix < numTSteps; ix++){
////            if ((LatLonAltStorage[px][ix][2] < NAS) && (shortestTime > ix)){
////                shortestTime = ix;
////                cout << "shortestTime = " << shortestTime << endl;
////                
////            }
////            //            cout << "Alt = " << LatLonAltStorage[0][ix][2] << endl;
////            
////            
////        }
////    }
//    
//    
//    for (int runx = 0; runx < num_per_batch; runx++){
//        int numTSteps = RocketLatLonAltStorage[runx].size();
//        cout << "runx = " << runx << endl;
//        
//        for (int ix = 0; ix < numTSteps; ix++){
//            if ((RocketLatLonAltStorage[runx][ix][2] < cutoff) && (longestTime < ix)){
//                longestTime = ix;
//                cout << "cur alt = " << RocketLatLonAltStorage[runx][ix][2] << endl;
//                cout << "longestTimeStep = " << longestTime << endl;   }
//        }
//    }
//    
//    cout << "longestTime in minutes = " << RocketStateTimes[longestTime]/60. << endl;
//    
//    timestepsUntilCutoff = longestTime;
//
//
//	return;
//}
//
//int Trajectory::getTimestepsUntilCutoff(){
//    return timestepsUntilCutoff;
//}
//
//
//
//
//
//void Trajectory::Propagate_First_Stage() {
//	
//	if (Properly_Initialized == false) {
//		cerr << "FATAL ERROR: You are trying to propagate but have not initialized yet!" << endl;
//		exit(666); }
//	
//	//Rocket Initial Conditions: Stationary on the launch pad
//	ColVec a_R_s__ecef(3,0.);	//Position of rocket (s) relative to launch pad (a) is zero vector
//	ColVec ecef_V_s__ecef(3,0.);	//Velocity of rocket (s) relative to launch pad (a) is zero vector
//	
//	ColVec NoUnderscores(3,0.);	//TEST
//
//	// Cape Canaveral coords
//	//WGS84	28° 23′ 18″ N, 80° 36′ 13″ W
//	double gdlat = 28.445455;	// geodetic degrees latitude
//	double lon = -80.564865;	// degrees longitude
//	double alt = 3e-3;			// altitude km above sea level
//	
//	ColVec o_R_a__ecef(3);
//	ColVec eci_omega_ecef__eci(3);
//
//	o_R_a__ecef = Geodetic_To_ECEF(gdlat, lon, alt);
//	//	cout << o_R_a__ecef;
//
//	eci_omega_ecef__eci(0) = 0;
//	eci_omega_ecef__eci(1) = 0;	
//	eci_omega_ecef__eci(2) = rotEarthRad;	
//	//	cout << eci_omega_ecef__eci;
//		
//	double UTC = Initial_UTC;
//	MatrixM eci__C__ecef(3,3);
//	eci__C__ecef = ECEF_to_ECI_rotation(UTC);
//	
//	
//	//Initial Conditions in ECI
//	ColVec o_R_s__eci(3,0.);	
//	ColVec eci_V_s__eci(3,0.);
//	
//	o_R_s__eci = prod(eci__C__ecef,(o_R_a__ecef + a_R_s__ecef));
//	eci_V_s__eci = prod(eci__C__ecef,ecef_V_s__ecef) + Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci);
//
////	cout << "CHEATING!  Using typed in values for initial position" << endl;
////	double temp1[3] = { -2.935341620790401e3,  -4.774747152374418e3,   3.029060975574545e3}; //[km]
////	o_R_s__eci(0) = temp1[0];
////	o_R_s__eci(1) = temp1[1];
////	o_R_s__eci(2) = temp1[2];
//
////	double temp2[3] = {0.3481800938695975,  -0.2140485115652455,   0};
////	eci_V_s__eci = temp2;
//	
//	initial_radius = VecNorm(o_R_s__eci);
//	cout << "initial radius = " << std::setprecision(16) << initial_radius << endl;
//
//	ColVec LatLonAlt(3,0.);
//	LatLonAlt = ECEF_To_Geodetic( o_R_s__eci);
//	
//	cout << "o_R_a__ecef = \n" <<  o_R_a__ecef;
//	cout << "o_R_s__eci = \n" <<  o_R_s__eci;
//	cout << "eci_V_s__eci = \n" << eci_V_s__eci;
//	cout << "Thetag = " << Sidereal_Time(Initial_UTC) << endl;
//	cout << "LatLonAlt = " << LatLonAlt[1]*180/PI << endl;
//
//	
//	
//	ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
//	state0(0) = o_R_s__eci(0);
//	state0(1) = o_R_s__eci(1);
//	state0(2) = o_R_s__eci(2);
//	state0(3) = eci_V_s__eci(0);
//	state0(4) = eci_V_s__eci(1);
//	state0(5) = eci_V_s__eci(2);
//	state0(6) = Mass0;
//	
//	cout << "state0 = " << state0 << endl;
//	cout << "Initial UTC = " << UTC << endl;
//
//	// Let's just make sure that we've got an initial state vector which makes sense
//	//   Convert eci_V_s__eci into local velocity
//	cout << "lat = " << gdlat << endl;
//	cout << "lon = " << lon << endl;
//	cout << "alt = " << alt << endl;
//	cout << "initial radius = " << std::setprecision(16) << initial_radius << endl;
//
//	// TODO: Why is this here?!  Already calculated this I think.
//	ColVec o_R_s__ecef(3,0.);
//	o_R_s__ecef = Geodetic_To_ECEF(gdlat, lon, alt);
////	cout << "o_R_s__ecef = " << o_R_s__ecef << endl;
//
//	MatrixM Local_C_ECEF (3,3,0.);
//	Local_C_ECEF = prod(Rotation(2, -gdlat), Rotation(3, lon));
////	cout << "Local_C_ECEF = " << Local_C_ECEF << endl;
//
//	MatrixM ECI__C__ECEF (3,3,0.);
//	ECI__C__ECEF = ECEF_to_ECI_rotation(UTC);
////	cout << "ECI__C__ECEF = " << ECI__C__ECEF << endl;
////	cout << "Initial UTC = " << UTC << endl;
//
//
//
////	ColVec o_R_s__local(3,0.);
////	o_R_s__local = prod(Local_C_ECEF,o_R_s__ecef);
//
////	cout << "o_R_s__eci = " << o_R_s__eci << endl;
//
//	ColVec eci_omega_local(3,0.);
//	eci_omega_local(2) = rotEarthRad;
////	cout << "eci_omega_local = " << eci_omega_local << endl;
//
//
//
//	MatrixM ECEF__C__Local(3,3,0.);
//	ECEF__C__Local = trans(Local_C_ECEF);
////	cout << "ECEF__C__Local = " << ECEF__C__Local << endl;
//
//	MatrixM ECI__C__Local;
//	ECI__C__Local = prod(ECI__C__ECEF, ECEF__C__Local);
////	cout << "ECI__C__Local = " << ECI__C__Local << endl;
//
//	ColVec Vlocal(3,1.);
//	Vlocal = prod(trans( ECI__C__Local),   eci_V_s__eci - Cross_Vectors(eci_omega_local, o_R_s__eci));
////	cout << "Vlocal = " << Vlocal << endl;
//
//
////	ColVec Vecef(3,0.);
////	Vecef = prod(ECEF__C__ECI,eci_V_s__eci);
////
////
////
////	ColVec Vlocal(3,0.);
////	Vlocal = prod(Local_C_ECEF,Vecef);
//
//	//-------------- End Debugging Stuff
//
//
//	
//	Initialize_density_profiles(num_per_batch);
//
////	timer stopwatch;
//	
////	stopwatch.start();
//
//	for (int i = 0; i < num_per_batch; i++) {
//		
//		load_density_profile(i);
//
//		gsl_integrate_first_stage(state0, i);
//		
//	}
////	stopwatch.stop();
//
//	
////	cout << "that took " << difftime(end_time,start_time) << " seconds" << endl;
////	cout << "in propagate(), that took " << stopwatch.how_long() << " seconds" << endl;
//	
//	
////	write_batch_to_file();
//    
////	write_batch_to_google_earth("meaningless",5);
//	
//	return;
//}
//
//void Trajectory::Propagate_Suborbital() {
//	
//	if (Properly_Initialized == false) {
//		cerr << "FATAL ERROR: You are trying to propagate but have not initialized yet!" << endl;
//		exit(666); }
//	
//    //Initial velocity of rocket in local frame
//    ColVec local_V_s__local(3,0.);
//    local_V_s__local(0) = Vx0*1e-3; //Convert to km/s
//    local_V_s__local(1) = Vy0*1e-3;
//    local_V_s__local(2) = Vz0*1e-3;
//    
//    //Rotate to ECEF (which is fixed wrt to local)
//    MatrixM Local_C_ECEF (3,3,0.);
//	Local_C_ECEF = prod(Rotation(2, -gdlat * PI/180), Rotation(3, lon * PI/180));
//    
//    MatrixM ECEF__C__Local(3,3,0.);
//	ECEF__C__Local = trans(Local_C_ECEF);
//    //	cout << "ECEF__C__Local = " << ECEF__C__Local << endl;
//    
//    
//	//Rocket Initial Conditions: Position is zero because launched from pad / plane
//	ColVec a_R_s__ecef(3,0.);	//Position of rocket (s) relative to launch pad (a) is zero vector
//        //Stays set to zero because rocket launches from pad (or plane)
//    
//    //Rocket Initial Conditions: Initial velocity of rocket due to being dropped from plane
//	ColVec ecef_V_s__ecef(3,0.);	//Velocity of rocket (s) relative to launch pad (a)
//	ecef_V_s__ecef = prod(ECEF__C__Local, local_V_s__local);  //(note: reference frame local fixed to ECEF)
//    
////    cout << "atan2() = " << atan2(local_V_s__local(2), local_V_s__local(1)) << endl;
//    
//    
//	// Distance from origin to launch pad / plane
//    local_elevation = 1.4;   //[km]
//    
//	ColVec o_R_a__ecef(3);  
//    o_R_a__ecef = Geodetic_To_ECEF(gdlat, lon, alt + local_elevation);
//
//    //Rotation of the earth
//	ColVec eci_omega_ecef__eci(3);
//	eci_omega_ecef__eci(0) = 0;
//	eci_omega_ecef__eci(1) = 0;
//	eci_omega_ecef__eci(2) = rotEarthRad;
//    
//    //Rotation matrix for converting from ecef to eci
//	double UTC = Initial_UTC;
//	MatrixM eci__C__ecef(3,3);
//	eci__C__ecef = ECEF_to_ECI_rotation(UTC);
//	
//	
//	//Initial Conditions in ECI
//	ColVec o_R_s__eci(3,0.);
//    o_R_s__eci = prod(eci__C__ecef,(o_R_a__ecef + a_R_s__ecef));
//    
//	ColVec eci_V_s__eci(3,0.);
//	eci_V_s__eci = prod(eci__C__ecef,ecef_V_s__ecef) + Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci);
//    
//    //~~~~~~~~~ Debug
//    cout << "eci_V_s__eci = " << eci_V_s__eci << endl;
//    cout << "ecef_V_s__eci = " << prod(eci__C__ecef,ecef_V_s__ecef) << endl;
//    cout << "CrossVec = " << Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci) << endl;
//    cout << "ecef_V_s__ecef = " << ecef_V_s__ecef << endl;
//    cout << "local_V_s__local = " << local_V_s__local << endl;
//    
//    
//    //Initial radius of earth gets used to calculate if we've crashed or not
//    //  This is lacking any elevation information.  Pretends earth is oblate spheroid.
//    initial_radius = get_local_earth_radius(gdlat*PI/180) + local_elevation;
//    
//    
//    //I'm guessing this was just to check that the conversions were working as expected
////	ColVec LatLonAlt(3,0.);
////	LatLonAlt = ECEF_To_Geodetic( o_R_s__eci);
////	
////	cout << "o_R_a__ecef = \n" <<  o_R_a__ecef;
////	cout << "o_R_s__eci = \n" <<  o_R_s__eci;
////	cout << "eci_V_s__eci = \n" << eci_V_s__eci;
////	cout << "Thetag = " << Sidereal_Time(Initial_UTC) << endl;
////	cout << "LatLonAlt = " << LatLonAlt[1]*180/PI << endl;
//    
//    cout << "~~~~~~~~~~~~~~~~~ Launch Info ~~~~~~~~~~~~~~~~~~~" << endl;
//    cout << "UTC = " << UTC << "   siderealTime = " << Sidereal_Time(Initial_UTC) << "  initialRadius = " << initial_radius << endl;
//    cout << "gdLat = " << gdlat << "   Lon = " << lon << "  Alt = " << alt << endl;
//    cout << "Vx0 = " << Vx0 << "  Vy0 = " << Vy0 << "  Vz0 = " << Vz0 << endl;
//    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
//    
//	
//	// Load up the initial state vector
//	ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
//	state0(0) = o_R_s__eci(0);
//	state0(1) = o_R_s__eci(1);
//	state0(2) = o_R_s__eci(2);
//	state0(3) = eci_V_s__eci(0);
//	state0(4) = eci_V_s__eci(1);
//	state0(5) = eci_V_s__eci(2);
//	state0(6) = Mass0;
//	
////	cout << "state0 = " << state0 << endl;
////	cout << "Initial UTC = " << UTC << endl;
////    
////	// Let's just make sure that we've got an initial state vector which makes sense
////	//   Convert eci_V_s__eci into local velocity
////	cout << "lat = " << gdlat << endl;
////	cout << "lon = " << lon << endl;
////	cout << "alt = " << alt << endl;
////	cout << "initial radius = " << std::setprecision(16) << initial_radius << endl;
//    
//	// TODO: Why is this here?!  Already calculated this I think.
////	ColVec o_R_s__ecef(3,0.);
////	o_R_s__ecef = Geodetic_To_ECEF(gdlat, lon, alt);
////	cout << "o_R_s__ecef = " << o_R_s__ecef << endl;
////
////	MatrixM ECI__C__ECEF (3,3,0.);
////	ECI__C__ECEF = ECEF_to_ECI_rotation(UTC);
////	cout << "ECI__C__ECEF = " << ECI__C__ECEF << endl;
////	cout << "Initial UTC = " << UTC << endl;
////    
////	cout << "o_R_s__eci = " << o_R_s__eci << endl;
////    
////	ColVec eci_omega_local(3,0.);
////	eci_omega_local(2) = rotEarthRad;
////	cout << "eci_omega_local = " << eci_omega_local << endl;
////    
////    
////	MatrixM ECI__C__Local;
////	ECI__C__Local = prod(ECI__C__ECEF, ECEF__C__Local);
////	cout << "ECI__C__Local = " << ECI__C__Local << endl;
////    
////	ColVec Vlocal(3,1.);
////	Vlocal = prod(trans( ECI__C__Local),   eci_V_s__eci - Cross_Vectors(eci_omega_local, o_R_s__eci));
////	cout << "Vlocal = " << Vlocal << endl;
//    
//	
//	Initialize_density_profiles(num_per_batch);
//    
//    //	timer stopwatch;
//	
//    //	stopwatch.start();
//    
//    int stateVecMaxSteps = 0;
//	for (int i = 0; i < num_per_batch; i++) {
//		
//		load_density_profile(i);
//        
//		int IX = gsl_integrate_suborbital(state0, i);
//        
////        Debris_State_Vector_Storage_Vec[piece_number].assign(IX,ZeroVecStateSize);
////		for (int i = 0; i < IX; i++) {
////			Debris_State_Vector_Storage_Vec[piece_number][i] = State_Vector_Buffer_Vec[i]; }
////        
//        
////        // Store the resulting state vectors
////        FallingStageStateVectorStorage[runNumber] = State_Vector_Buffer_Vec;
////        FallingStageStateVectorStorage[runNumber].resize(IX);
////        
////        //        State_Vector_Buffer_Vec.resize(IX);
////        //        StageDownStateTimes.resize(IX);
////        
//        // debris_max_time_steps keeps track of the maximum number of time steps used
//		if (IX > stateVecMaxSteps) {
//			stateVecMaxSteps = IX; }
////    }
////    
//		
//        
//        
//	}
//    //	stopwatch.stop();
//    
//    // Truncate state times to appropriate length
//    SuborbitalStateTimes.resize(stateVecMaxSteps);
//
//    
//	
//    //	cout << "that took " << difftime(end_time,start_time) << " seconds" << endl;
//    //	cout << "in propagate(), that took " << stopwatch.how_long() << " seconds" << endl;
//	
//	
////	write_batch_to_file();
////	write_batch_to_google_earth("meaningless",5);
//	
//	return;
//}
//
//
//
//
//
//void Trajectory::Initialize_density_profiles(int num_profiles) {
////	cout << "num_profiles = " << num_profiles << endl;
//	
//	switch (Density_Option) {
//		case VACUUM:
//		case GRAM_DENSITY:
//		case CANTWELL_DENSITY:
//			// Do nothing and like it
//			break;
//			
//		case GRAM_DENSITY_UNCERTAIN: {			
//			double delta_step = (max_alt_rho - min_alt_rho)/(num_uncert_density_steps - 1.);
//
//			// Allocate density profile altitude array
//			Density_uncert_steps = new double[num_uncert_density_steps];
//			for (int i = 0; i < num_uncert_density_steps; i++) {
//				Density_uncert_steps[i] = min_alt_rho + i*delta_step; }
//			
//			// Allocate density profiles matrix
//			Density_profiles = new double* [num_profiles];
//			for (int i = 0; i < num_profiles; i++) {
//				Density_profiles[i] = new double[num_uncert_density_steps]; }
//			
//			// Initialize density profiles
//			for (int alt_j = 0; alt_j < num_uncert_density_steps; alt_j++) {
//				double rho_mean = gsl_spline_eval (RhoMean_spline, Density_uncert_steps[alt_j], RhoMean_acc);
//				double rho_std_percent = gsl_spline_eval (RhoSd_spline, Density_uncert_steps[alt_j], RhoSd_acc);
//				
//				for (int i = 0; i < num_profiles; i++) {
//					
//					double random_density = abs( rand_num->generate_random_gaussian(rho_mean, rho_std_percent*rho_mean/100) );
//					//cout << "random_density = " << abs(rand_num->generate_random_gaussian(5,1)) << endl;
//					Density_profiles[i][alt_j] = random_density;	} }
//			} break;
//	
//		
//		default:
//			cout << "ERROR!!! Initialize density profiles went wrong!!!~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//			break;
//	}
//	
//	return;
//}
//
//
//
//void Trajectory::load_density_profile(int run_number) {
//	
//	
//	switch (Density_Option) {
//		case GRAM_DENSITY_UNCERTAIN:
////			cout << "here with i = " << i << endl;
//			
//			//Set up the interpolators (GSL cubic spline)
//			RhoUncert_acc = gsl_interp_accel_alloc ();
//			RhoUncert_spline = gsl_spline_alloc (gsl_interp_linear, num_uncert_density_steps);
//			gsl_spline_init (RhoUncert_spline, Density_uncert_steps, Density_profiles[run_number], num_uncert_density_steps);
//			
//			
//			break;
//			
//		case CANTWELL_DENSITY:
//		case GRAM_DENSITY:
//		case VACUUM:
//			//Do nothing and like it
//			break;
//			
//		default:
//			cout << "ERROR in load density profile! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//		break; }
//	
//
//	
//	return;
//}


//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Vehicle Properties and Setup Files
//
//void Trajectory::Get_Nominal_Trajectory(string nominal_file) {
//	
//	ifstream opt_traj;
//	
//	//read in the optimal trajectory file.  MUST USE c_str() or else this will fail!!!
//	opt_traj.open(nominal_file.c_str(), ifstream::in);
//	
//	if (opt_traj.is_open()){
////		cout << "optimal trajectory file opened successfully!\n";	
//	}
//	else {
//		std::cerr << "optimal trajectory file failed to open :(\n";	}
//	
//	int buf_size = 500;
//	char buffer[buf_size];
//	//Want to kill one line
//	for (int i=0; i < 1; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	opt_traj >> Fmax >> Isp >> Mass0 >> NumSteps >> NumStages;
//	
////	cout << "initial mass = " << Mass0 << endl;
////	cout << "NumStages = " << NumSteps << endl;
//	
//	StagePropMass = new double[NumStages];
//	StageStructMass = new double[NumStages];
//	StageTimes = new double[NumStages];
//	StageISP = new double[NumStages];
//	StageFmax = new double[NumStages];
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageStructMass[i]; }
//	
////	cout << "StagePropMass 2 = " << StageStructMass[1] << endl;
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StagePropMass[i]; }
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageISP[i]; }	
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageFmax[i]; }		
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageTimes[i]; }	
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	opt_traj >> PayloadMass;
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	opt_traj >> gdlat >> lon >> alt;
//	
////	cout << "initial longitude = " << lon << endl;
//	
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	opt_traj >> launch_day >> launch_month >> launch_year >> launch_hours >> launch_minutes >> timezone_offset;
//	
////	cout << "launch year = " << launch_year << endl;
//	
//	Initial_UTC = UTC_Time(launch_day, launch_month, launch_year, launch_hours, launch_minutes, timezone_offset);
//	cout << "launch_UTC = " << Initial_UTC << endl;
//	
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//	
//	//allocate
//	ThrustUx = new double[NumSteps];
//	ThrustUy = new double[NumSteps];
//	ThrustUz = new double[NumSteps];
//
//	ThrustMag = new double[NumSteps];
//	ThrustTime = new double[NumSteps];
//	
//	for (int i=0; i < NumSteps; i++) {
//		opt_traj >> ThrustUx[i] >> ThrustUy[i] >> ThrustUz[i] >> ThrustMag[i] >> ThrustTime[i]; }
//	
//	for (int i=0; i < NumSteps; i++) {
//		cout << ThrustUx[i] << "  " << ThrustUy[i] << "  " << ThrustUz[i] << "  " << ThrustMag[i] << "  " << ThrustTime[i] << endl; }
//
//	opt_traj.close();
//	
//	//save thrust time limits
//	min_thrust_time = ThrustTime[0];
//	max_thrust_time = ThrustTime[NumSteps-1];
//	
//	
//	//NOTE: Splines have problems when Mag goes to zero because you get degenerate time values.  Must be fixed for staging!!!
//	//Set up the interpolators (GSL cubic spline)
//	ThrustUx_acc = gsl_interp_accel_alloc ();
//	ThrustUx_spline = gsl_spline_alloc (gsl_interp_linear, NumSteps);
//	gsl_spline_init (ThrustUx_spline, ThrustTime, ThrustUx, NumSteps);
//	
//	//Set up the interpolators (GSL cubic spline)
//	ThrustUy_acc = gsl_interp_accel_alloc ();
//	ThrustUy_spline = gsl_spline_alloc (gsl_interp_linear, NumSteps);
//	gsl_spline_init (ThrustUy_spline, ThrustTime, ThrustUy, NumSteps);
//	
//	//Set up the interpolators (GSL cubic spline)
//	ThrustUz_acc = gsl_interp_accel_alloc ();
//	ThrustUz_spline = gsl_spline_alloc (gsl_interp_linear, NumSteps);
//	gsl_spline_init (ThrustUz_spline, ThrustTime, ThrustUz, NumSteps);
//	
//	//Set up the interpolators (GSL cubic spline)
//	ThrustMag_acc = gsl_interp_accel_alloc ();
//	ThrustMag_spline = gsl_spline_alloc (gsl_interp_linear, NumSteps);
//	gsl_spline_init (ThrustMag_spline, ThrustTime, ThrustMag, NumSteps);	
//	
//	
//	//Destructor sometime...right now this is a big memory leak
//	
//	return;
//}
//
//
//void Trajectory::Get_Nominal_Suborbital(string nominal_file) {
//
//	ifstream opt_traj;
//
//	//read in the optimal trajectory file.  MUST USE c_str() or else this will fail!!!
//	opt_traj.open(nominal_file.c_str(), ifstream::in);
//
//	if (opt_traj.is_open()){
////		cout << "optimal trajectory file opened successfully!\n";
//	}
//	else {
//		std::cerr << "optimal trajectory file failed to open :(\n";	}
//
//	int buf_size = 500;
//	char buffer[buf_size];
//	//Want to kill one line
//	for (int i=0; i < 1; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	opt_traj >> Fmax >> Isp >> Mass0 >> NumSteps >> NumStages;
//
////	cout << "initial mass = " << Mass0 << endl;
////	cout << "NumStages = " << NumSteps << endl;
//
//	StagePropMass = new double[NumStages];
//	StageStructMass = new double[NumStages];
//	StageTimes = new double[NumStages];
//	StageISP = new double[NumStages];
//	StageFmax = new double[NumStages];
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageStructMass[i]; }
//
////	cout << "StagePropMass 2 = " << StageStructMass[1] << endl;
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StagePropMass[i]; }
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageISP[i]; }
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageFmax[i]; }
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < NumStages; i++) {
//		opt_traj >> StageTimes[i]; }
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	opt_traj >> PayloadMass;
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	opt_traj >> gdlat >> lon >> alt;
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	// Load initial velocities in local frame
//	opt_traj >> Vx0 >> Vy0 >> Vz0;
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	opt_traj >> launch_day >> launch_month >> launch_year >> launch_hours >> launch_minutes >> timezone_offset;
//
////	cout << "launch year = " << launch_year << endl;
//
//	Initial_UTC = UTC_Time(launch_day, launch_month, launch_year, launch_hours, launch_minutes, timezone_offset);
//	cout << "launch_UTC = " << Initial_UTC << endl;
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	// Load initial velocities in local frame
//	opt_traj >> numAccX >> numAccY >> numAccZ;
//
//	//allocate
//	double *accNX = new double[numAccX];
//	double *accNY = new double[numAccY];
//	double *accNZ = new double[numAccZ];
//
//	double *timeNX = new double[numAccX];
//	double *timeNY = new double[numAccY];
//	double *timeNZ = new double[numAccZ];
//
////	ThrustMag = new double[NumSteps];
////	ThrustTime = new double[NumSteps];
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < numAccX; i++) {
//		opt_traj >> timeNX[i] >> accNX[i]; }
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < numAccY; i++) {
//		opt_traj >> timeNY[i] >> accNY[i]; }
//
//	//Want to finish this line and kill two more
//	for (int i=0; i < 3; i++){
//		opt_traj.getline(buffer,buf_size);	}
//
//	for (int i=0; i < numAccZ; i++) {
//		opt_traj >> timeNZ[i] >> accNZ[i]; }
//
////	for (int i=0; i < NumSteps; i++) {
////		cout << ThrustUx[i] << "  " << ThrustUy[i] << "  " << ThrustUz[i] << "  " << ThrustMag[i] << "  " << ThrustTime[i] << endl; }
//
//	opt_traj.close();
//
//	//save thrust time limits
//    NXstart = timeNX[0];
//    NXend = timeNX[numAccX-1];
//    NYstart = timeNY[0];
//    NYend = timeNY[numAccY-1];
//    NZstart = timeNZ[0];
//    NZend = timeNZ[numAccZ-1];
//
//
//	//NOTE: Splines have problems when Mag goes to zero because you get degenerate time values.  Must be fixed for staging!!!
//	//Set up the interpolators (GSL cubic spline)
//	ThrustUx_acc = gsl_interp_accel_alloc ();
//	ThrustUx_spline = gsl_spline_alloc (gsl_interp_linear, numAccX);
//	gsl_spline_init (ThrustUx_spline, timeNX, accNX, numAccX);
//
//	//Set up the interpolators (GSL cubic spline)
//	ThrustUy_acc = gsl_interp_accel_alloc ();
//	ThrustUy_spline = gsl_spline_alloc (gsl_interp_linear, numAccY);
//	gsl_spline_init (ThrustUy_spline, timeNY, accNY, numAccY);
//
//	//Set up the interpolators (GSL cubic spline)
//	ThrustUz_acc = gsl_interp_accel_alloc ();
//	ThrustUz_spline = gsl_spline_alloc (gsl_interp_linear, numAccZ);
//	gsl_spline_init (ThrustUz_spline, timeNZ, accNZ, numAccZ);
//
//	//Destructor sometime...right now this is a big memory leak
//    
//    
//    //Calc initial heading angle wrt to local y (east)
//    headingAngle0 = atan2(Vz0, Vx0);
//    
//    
//    
//    
//    
//    
//    
//    //While we're at it, let's also read in the pitching model
//    string pitching_file = "Files/SS2/SS2PitchModel.txt";
//    ifstream pitchModel;
//    
//	//read in the optimal trajectory file.  MUST USE c_str() or else this will fail!!!
//	pitchModel.open(pitching_file.c_str(), ifstream::in);
//    
//	if (pitchModel.is_open()){
//        //		cout << "optimal trajectory file opened successfully!\n";
//	}
//	else {
//		std::cerr << "pitch model file failed to open :(\n";	}
//    
//    //Want to kill three lines
//	for (int i=0; i < 3; i++){
//		pitchModel.getline(buffer,buf_size); }
//
//    int numStepsPitch;
//    pitchModel >> numStepsPitch;
////    numStepsPitch = 60;
//    
//    //Finish this line and to kill one more
//	for (int i=0; i < 2; i++){
//		pitchModel.getline(buffer,buf_size); }
//
//    double *PitchTime = new double[numStepsPitch];
//	double *PitchValue = new double[numStepsPitch];
//	double *OmegaValue = new double[numStepsPitch];
//	double *OmegaDotValue = new double[numStepsPitch];
//    
//    for (int i=0; i < numStepsPitch; i++) {
//        pitchModel >> PitchTime[i] >> PitchValue[i] >> OmegaValue[i] >> OmegaDotValue[i]; }
//
//    pitchModel.close();
//    
//    for (int i = 0; i < numStepsPitch; i++){
//        cout << "X = " << PitchTime[i] << endl;
//    }
//    
//    
//    //Set up the interpolators (GSL cubic spline)
//	PitchAngle_acc = gsl_interp_accel_alloc ();
//	PitchAngle_spline = gsl_spline_alloc (gsl_interp_linear, numStepsPitch);
//    Omega_acc = gsl_interp_accel_alloc ();
//	Omega_spline = gsl_spline_alloc (gsl_interp_linear, numStepsPitch);
//    OmegaDot_acc = gsl_interp_accel_alloc ();
//	OmegaDot_spline = gsl_spline_alloc (gsl_interp_linear, numStepsPitch);
//    
//	gsl_spline_init (PitchAngle_spline, PitchTime, PitchValue, numStepsPitch);
//	gsl_spline_init (Omega_spline, PitchTime, OmegaValue, numStepsPitch);
//	gsl_spline_init (OmegaDot_spline, PitchTime, OmegaDotValue, numStepsPitch);
//    
//    //Save the time interval
//    timeAngleStart = PitchTime[0];
//    timeAngleEnd = PitchTime[numStepsPitch-1];
//    
//    
//    
//    
//
//	return;
//}
//
//
//
//
//
//void Trajectory::Get_Approximate_Wind(string approx_wind_filename) {
//	
//	ifstream wind_file;
//	
//	//read in the optimal trajectory file.  MUST USE c_str() or else this will fail!!!
//	wind_file.open(approx_wind_filename.c_str(), ifstream::in);
//	
//	if (wind_file.is_open()){
//				cout << "wind file file opened successfully!\n";	
//	}
//	else {
//		std::cerr << "GRAM wind file failed to open :(\n File was called: " << approx_wind_filename << endl;	}
//	
//	char buffer[500];
//	//Want to kill one line
//	for (int i=0; i < 1; i++){
//		wind_file.getline(buffer,500);}
//	
//	int NumPoints;
//	wind_file >> NumPoints;
//	
//	//Want to finish current line and kill two more lines
//	for (int i=0; i < 3; i++){
//		wind_file.getline(buffer,500);}
//	
//	cout << "Number of data points = " << NumPoints << endl;
//	
//	double *RhoMean = new double[NumPoints];
//	double *WindUMean = new double[NumPoints];
//	double *WindVMean = new double[NumPoints];
//	double *WindWMean = new double[NumPoints];
//	double *RhoSd = new double[NumPoints];
//	double *WindUSd = new double[NumPoints];
//	double *WindVSd = new double[NumPoints];
//	double *WindWSd = new double[NumPoints];
//	
//	double *WindAlt = new double[NumPoints];
//	
//	cout << "THE WAY THE WIND FILE IS SPECIFIED SUCKS, FIX THIS LATER" << endl;
//	
//	//note, i'm loading these arrays backwards so that the altitude array will be monotonically increasing
//	for (int i = NumPoints-1; i >= 0; i--) {
//		wind_file >> WindAlt[i] >> RhoMean[i] 
//		>> WindUMean[i] >> WindVMean[i] >> WindWMean[i] >> RhoSd[i] >> WindUSd[i] >> WindVSd[i] >> WindWSd[i]; 
//		
////		cout << i << "  wind alt = " << WindAlt[i] << endl;
//	}
//	
//	wind_file.close();
//	
//	//store altitude limits
//	min_alt_rho = WindAlt[0];
//	max_alt_rho = WindAlt[NumPoints-1];
//	min_alt_wind = WindAlt[0];
//	max_alt_wind = WindAlt[NumPoints-1];
//	
//	//Set up the interpolators (GSL cubic spline)
//	RhoMean_acc = gsl_interp_accel_alloc ();
//	RhoMean_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (RhoMean_spline, WindAlt, RhoMean, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	WindUMean_acc = gsl_interp_accel_alloc ();
//	WindUMean_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (WindUMean_spline, WindAlt, WindUMean, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	WindVMean_acc = gsl_interp_accel_alloc ();
//	WindVMean_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (WindVMean_spline, WindAlt, WindVMean, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	WindWMean_acc = gsl_interp_accel_alloc ();
//	WindWMean_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (WindWMean_spline, WindAlt, WindWMean, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	RhoSd_acc = gsl_interp_accel_alloc ();
//	RhoSd_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (RhoSd_spline, WindAlt, RhoSd, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	WindUSd_acc = gsl_interp_accel_alloc ();
//	WindUSd_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (WindUSd_spline, WindAlt, WindUSd, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	WindVSd_acc = gsl_interp_accel_alloc ();
//	WindVSd_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (WindVSd_spline, WindAlt, WindVSd, NumPoints);
//	
//	//Set up the interpolators (GSL cubic spline)
//	WindWSd_acc = gsl_interp_accel_alloc ();
//	WindWSd_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (WindWSd_spline, WindAlt, WindWSd, NumPoints);
//	
//
//	// Generate a vector of evenly spaced wind data for quick lookup and see if faster than interpolation
//		// Do this when you don't have much else better to do.
//	atm_bin_size_km = 1;
//	num_wind_bins = (INTxx) ceil(max_alt_wind/atm_bin_size_km);
//	WindProfileUVW.assign(num_wind_bins, ColVec(3,-5.));		//memory leak
//
//	initialize_random_atmosphere();
//
//	return;
//}
//
//int Trajectory::initialize_random_atmosphere() {
//	//Note that wind units are m/s!!!  Alt units are km!!!
//	//want to do a random walk to get things continuous
//
//	Random_Number myrand;
//
//	double max_u_slope = 20; // set to about twice the maximum possible slope (m/s)/km
//	double max_v_slope = 12; //
//	double max_w_slope = 3; //
//	// Just for fun, let's generate a bunch of distributions and plot them against the actual envelopes
//	int num_to_generate = 1;
//
//	for (int jx = 0; jx < num_to_generate; jx++) {
//
//		//Pick first point to get things started
//		double alt = 0;
//		double u_mean = gsl_spline_eval (WindUMean_spline, alt, WindUMean_acc); //East
//		double v_mean = gsl_spline_eval (WindVMean_spline, alt, WindVMean_acc); //North
//		double w_mean = gsl_spline_eval (WindWMean_spline, alt, WindWMean_acc); //Up
//
//		double u_std = gsl_spline_eval (WindUSd_spline, alt, WindUSd_acc);
//		double v_std = gsl_spline_eval (WindVSd_spline, alt, WindVSd_acc);
//		double w_std = gsl_spline_eval (WindWSd_spline, alt, WindWSd_acc);
//
//		WindProfileUVW[0][0] = myrand.generate_random_gaussian(u_mean, u_std);
//		WindProfileUVW[0][1] = myrand.generate_random_gaussian(v_mean, v_std);
//		WindProfileUVW[0][2] = myrand.generate_random_gaussian(w_mean, w_std);
//
//		// Step size will be std dev of previous point scaled by size of altitude bins
//		double u_step_size = max_u_slope*atm_bin_size_km;
//		double v_step_size = max_v_slope*atm_bin_size_km;
//		double w_step_size = max_w_slope*atm_bin_size_km;
//
//		//Now random walk from there
//		for (int ix = 1; ix < num_wind_bins; ix++) {
//			alt = ix*atm_bin_size_km;
//			u_mean = gsl_spline_eval (WindUMean_spline, alt, WindUMean_acc); //East
//			v_mean = gsl_spline_eval (WindVMean_spline, alt, WindVMean_acc); //North
//			w_mean = gsl_spline_eval (WindWMean_spline, alt, WindWMean_acc); //Up
//
//			u_std = gsl_spline_eval (WindUSd_spline, alt, WindUSd_acc);
//			v_std = gsl_spline_eval (WindVSd_spline, alt, WindVSd_acc);
//			w_std = gsl_spline_eval (WindWSd_spline, alt, WindWSd_acc);
//
//			// These get used for all three directions U, V, and W
//			double chooseU = 1;
//			double chooseV = 1;
//			double chooseW = 1;
//			double prob_u_happens = 0;
//			double prob_v_happens = 0;
//			double prob_w_happens = 0;
//
//			// Look at u-direction
//			double u_step = -5;
//			double v_step = -5;
//			double w_step = -5;
//
//			while ( (chooseU > prob_u_happens) || (chooseV > prob_v_happens) || (chooseW > prob_w_happens)) {
//				// Find out what percent of the step_size me move and in what direction
//				double step_scale = myrand.generate_random_uniform(-1,1);
//
//				// take the proposed step
//				u_step = WindProfileUVW[ix-1][0] + step_scale*u_step_size;
//				v_step = WindProfileUVW[ix-1][1] + step_scale*v_step_size;
//				w_step = WindProfileUVW[ix-1][2] + step_scale*w_step_size;
//
//				// Check to see how likely the step is
////				double u_upper_bound = u_step - u_mean + u_step_size/30;
////				double u_lower_bound = u_step - u_mean - u_step_size/30;
////				double v_upper_bound = v_step - v_mean + v_step_size;
////				double v_lower_bound = v_step - v_mean - v_step_size;
////				double w_upper_bound = w_step - w_mean + w_step_size;
////				double w_lower_bound = w_step - w_mean - w_step_size;
//
////				prob_u_happens = myrand.gaussian_cdf(u_upper_bound, u_std) - myrand.gaussian_cdf(u_lower_bound, u_std);
//				prob_u_happens = myrand.gaussian_cdf(u_mean - u_step, u_std);
//				if (prob_u_happens > 0.5) {
//					prob_u_happens = 1 - prob_u_happens;
//				}
//
//				prob_v_happens = myrand.gaussian_cdf(v_mean - v_step, v_std);
//				if (prob_v_happens > 0.5) {
//					prob_v_happens = 1 - prob_v_happens;
//				}
//
//				prob_w_happens = myrand.gaussian_cdf(w_mean - w_step, w_std);
//				if (prob_w_happens > 0.5) {
//					prob_w_happens = 1 - prob_w_happens;
//				}
//
////				prob_v_happens = 1.;//myrand.gaussian_cdf(v_upper_bound, v_std) - myrand.gaussian_cdf(v_lower_bound, v_std);
////				prob_w_happens = 1.;//myrand.gaussian_cdf(w_upper_bound, w_std) - myrand.gaussian_cdf(w_lower_bound, w_std);
//
//				chooseU = myrand.generate_random_uniform(0.,1.);
//				chooseV = myrand.generate_random_uniform(0.,1.);
//				chooseW = myrand.generate_random_uniform(0.,1.);
////				cout << "choose = " << chooseU << " < " << prob_u_happens << endl;
//			}
////			cout << "found" << endl;
//
//			// Left the while loop so the current values are good
//			WindProfileUVW[ix][0] = u_step;
//			WindProfileUVW[ix][1] = v_step;
//			WindProfileUVW[ix][2] = w_step;
//
////			// Save std devs for next step
////			u_step_size = max_u_slope*atm_bin_size_km;
////			v_step_size = v_std*atm_bin_size_km;
////			w_step_size = w_std*atm_bin_size_km;
//		}
//
//		// If num_to_generate is greater than 1, i presume we're debugging or making plots, so print stuff out in that case
//		if (num_to_generate > 1) {
//			//Print the results to check in matlab
//			cout << "\n\n\n\n\nrandom_walk{" << jx+1 << "} = [\n";
//			for (int ix = 0; ix < num_wind_bins; ix++) {
//				alt = ix*atm_bin_size_km;
//				cout << "  " << alt << "   " << WindProfileUVW[ix][0] << "   " << WindProfileUVW[ix][1] << "   " << WindProfileUVW[ix][2] << endl;
//			}
//			cout << "];" << endl;
//		}
//	}
//
//
//	if (num_to_generate > 1){
//		cout << "alt_mean3_std3 = [\n";
//		for (int ix = 0; ix < num_wind_bins; ix++) {
//					alt = ix*atm_bin_size_km;
//					double u_mean = gsl_spline_eval (WindUMean_spline, alt, WindUMean_acc); //East
//					double v_mean = gsl_spline_eval (WindVMean_spline, alt, WindVMean_acc); //North
//					double w_mean = gsl_spline_eval (WindWMean_spline, alt, WindWMean_acc); //Up
//
//					double u_std = gsl_spline_eval (WindUSd_spline, alt, WindUSd_acc);
//					double v_std = gsl_spline_eval (WindVSd_spline, alt, WindVSd_acc);
//					double w_std = gsl_spline_eval (WindWSd_spline, alt, WindWSd_acc);
//
//					cout << "  " << alt << "  " << u_mean << "  " << v_mean << "  " << w_mean << "  " << u_std << "  " << v_std << "  " << w_std << endl;
//		}
//		cout << "];\n\n\n";
//		cout << "EXIT!!!!!!" << endl;
//		cout << "don't forget that you changed the V and W outputs to be an envelope around U" << endl;
//		exit(234);
//	}
//
//
//
//	return 1;
//}
//
//
//
//
//
//
//
//
//void Trajectory::Get_Temperature_Data(string temperature_file){
//	
//	ifstream temp_file;
//	
//	//read in the optimal trajectory file.  MUST USE c_str() or else this will fail!!!
//	temp_file.open(temperature_file.c_str(), ifstream::in);
//	
//	if (temp_file.is_open()){
//		cout << "temperature file opened successfully!\n";	
//	}
//	else {
//		std::cerr << "temperature file failed to open :(\n";	}
//	
//	char buffer[500];
//	//Want to kill one line
//	for (int i=0; i < 1; i++){
//		temp_file.getline(buffer,500);}
//	
//	int NumPoints;
//	temp_file >> NumPoints;
//	
//	//Want to finish current line and kill one more lines
//	for (int i=0; i < 2; i++){
//		temp_file.getline(buffer,500);}
//	
//	cout << "Number of data points = " << NumPoints << endl;
//	
//	double *TempProfile = new double[NumPoints];
//	double *AltTemp = new double[NumPoints];
//
//	
//	for (int i = 0; i < NumPoints; i++) {
//		temp_file >> AltTemp[i] >> TempProfile[i]; 
//	
////		cout << i << "  " << TempProfile[i] << endl; 
//	}
//	
//	temp_file.close();
//	
//	//store temp altitude limits
//	min_alt_temp = AltTemp[0];
//	max_alt_temp = AltTemp[NumPoints-1];
//	
//	//Set up the interpolators (GSL cubic spline)
//	Temp_acc = gsl_interp_accel_alloc ();
//	Temp_spline = gsl_spline_alloc (gsl_interp_linear, NumPoints);
//	gsl_spline_init (Temp_spline, AltTemp, TempProfile, NumPoints);
//	
//	return;
//}
//
//




















//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//ColVec Trajectory::GetThrustForce(double cur_time) {
//	ColVec ThrustForce(3,0.0);
//	
//	if ((cur_time >= min_thrust_time) && (cur_time < max_thrust_time)) {
//		double ux = gsl_spline_eval (ThrustUx_spline, cur_time, ThrustUx_acc);
//		double uy = gsl_spline_eval (ThrustUy_spline, cur_time, ThrustUy_acc);
//		double uz = gsl_spline_eval (ThrustUz_spline, cur_time, ThrustUz_acc);
//		double eta = gsl_spline_eval (ThrustMag_spline, cur_time, ThrustMag_acc);
//		
//		ThrustForce(0) = Fmax*eta*ux;	//currently in Newtons
//		ThrustForce(1) = Fmax*eta*uy;
//		ThrustForce(2) = Fmax*eta*uz;
//	}
//	else {
//		ThrustForce(0) = 0;
//		ThrustForce(1) = 0;
//		ThrustForce(2) = 0;
//	}
//	
//	return ThrustForce;
//}
//
//double Trajectory::SuborbitalTotalAcc(double cur_time){
//	
//    // X Direction
//    double xAcc = -5;
//	if ((cur_time >= NXstart) && (cur_time < NXend)) {
//        xAcc = gsl_spline_eval (ThrustUx_spline, cur_time, ThrustUx_acc); }
//	else {
//		xAcc = 0; }
//    
//    // For now I know that I'm not using this, so comment it out
//    //    // Y Direction
//    //	if ((cur_time >= NYstart) && (cur_time < NYend)) {
//    //        ThrustAcc(1) = gsl_spline_eval (ThrustUy_spline, cur_time, ThrustUy_acc); }
//    //	else {
//    //		ThrustAcc(1) = 0; }
//    
//    
//    // Z Direction - Per my convention, this will only be used to calculate the pitch of the vehicle
//    //   NO!  ALSO NEED FOR MDOT CALCULATION
//    double zAcc = -5;
//	if ((cur_time >= NZstart) && (cur_time < NZend)) {
//        zAcc = gsl_spline_eval (ThrustUz_spline, cur_time, ThrustUz_acc); }
//	else {
//		zAcc = 0; }
//    //Actually...i probably don't even need this.  just use the state time.
//    
//    double totalAcc = sqrt(xAcc*xAcc + zAcc*zAcc);
//	
//	return totalAcc;
//}
//
//
//ColVec Trajectory::GetBodyAcc(double cur_time) {
//	ColVec ThrustAcc(3,0.0);
//	
//    // X Direction
//    double xAcc = -5;
//	if ((cur_time >= NXstart) && (cur_time < NXend)) {
//        xAcc = gsl_spline_eval (ThrustUx_spline, cur_time, ThrustUx_acc); }
//	else {
//		xAcc = 0; }
//    ThrustAcc(0) = xAcc;
//    
//    // Y Direction - Ignoring for the moment
//    
//    // Z Direction
//    double zAcc = -5;
//	if ((cur_time >= NZstart) && (cur_time < NZend)) {
//        zAcc = gsl_spline_eval (ThrustUz_spline, cur_time, ThrustUz_acc); }
//	else {
//		zAcc = 0; }
//    ThrustAcc(2) = zAcc;
//    
//	return ThrustAcc;   // value is in g's
//}
//
//ColVec Trajectory::GetBodyAngles(double cur_time) {
//    ColVec bodyAngles(3,0.0);
//
//    double pitchAngle = -5;
//    if ((cur_time >= timeAngleStart) && (cur_time < timeAngleEnd)) {
//        pitchAngle = gsl_spline_eval (PitchAngle_spline, cur_time, PitchAngle_acc); }
//	else if (cur_time < timeAngleStart) {
//		pitchAngle = 0; }
//    else {
//        pitchAngle = gsl_spline_eval (PitchAngle_spline, timeAngleEnd, PitchAngle_acc); }
//    bodyAngles(0) = pitchAngle;  
//    
//    
//    double omega = -5;
//    if ((cur_time >= timeAngleStart) && (cur_time < timeAngleEnd)) {
//        omega = gsl_spline_eval (Omega_spline, cur_time, Omega_acc); }
//	else if (cur_time < timeAngleStart) {
//		omega = 0; }
//    else {
//        omega = gsl_spline_eval (Omega_spline, timeAngleEnd, Omega_acc); }
//    bodyAngles(1) = omega;
//    
//    
//    double OmegaDot = -5;
//    if ((cur_time >= timeAngleStart) && (cur_time < timeAngleEnd)) {
//        OmegaDot = gsl_spline_eval (OmegaDot_spline, cur_time, OmegaDot_acc); }
//	else if (cur_time < timeAngleStart) {
//		OmegaDot = 0; }
//    else {
//        OmegaDot = gsl_spline_eval (OmegaDot_spline, timeAngleEnd, OmegaDot_acc); }
//    bodyAngles(2) = OmegaDot;
//    
//    return bodyAngles;
//    
//}
//
//
//ColVec Trajectory::VwindCurveUVW(ColVec r__ECI) {
//	//U => east
//	//V => north
//	//W => up
//	
//	ColVec Vwind(3,0.);
//	
//	double altitude = VecNorm(r__ECI) - initial_radius;
//	
//	if ((altitude >= min_alt_wind) && (altitude < max_alt_wind)) {
//		double u = gsl_spline_eval (WindUMean_spline, altitude, WindUMean_acc); //East
//		double v = gsl_spline_eval (WindVMean_spline, altitude, WindVMean_acc); //North
//		double w = gsl_spline_eval (WindWMean_spline, altitude, WindWMean_acc); //Up
//		
//		Vwind(0) = u;	//East
//		Vwind(1) = -1000*v;	//North
//		Vwind(2) = w;	//Up
//
//	}
//	else if (altitude < min_alt_wind){
////		cout << "you may have crashed Vwind = 0" << endl; 
//	}
//	else {
//		//Vwind is zero vector
//	}
//	
//	cout << "min_alt_wind " << min_alt_wind << endl;
//	cout << "max_alt_wind " << max_alt_wind << endl;
//	cout << "initial_radius " << initial_radius << endl;
//	cout << "altitude " << altitude << endl;
//	cout << "Vwind = " << Vwind << endl;
//
//	return Vwind*1e-3;	//convert to km/s
//}
//
//ColVec Trajectory::VwindCurveUVW() {
//	// This version of the function doesn't take any inputs but rather draws from private class variables
//	//  This is kinda dangerous, but it'll work for now.
//
//	//U => east
//	//V => north
//	//W => up
//
//	ColVec Vwind(3,0.);
//
////	double altitude = current_alt;
//
//	int alt_index = (INTxx) floor(current_alt/atm_bin_size_km);
//	if ((alt_index < num_wind_bins) && (alt_index >= 0)) {
//		Vwind = WindProfileUVW[alt_index]; }
//
//
////	if ((altitude >= min_alt_wind) && (altitude < max_alt_wind)) {
////		double u = gsl_spline_eval (WindUMean_spline, altitude, WindUMean_acc); //East
////		double v = gsl_spline_eval (WindVMean_spline, altitude, WindVMean_acc); //North
////		double w = gsl_spline_eval (WindWMean_spline, altitude, WindWMean_acc); //Up
////
////		Vwind(0) = u;	//East
////		Vwind(1) = v;	//North
////		Vwind(2) = w;	//Up
////
////	}
////	else if (altitude < min_alt_wind){
//////		cout << "you may have crashed Vwind = 0" << endl;
////	}
////	else {
////		//Vwind is zero vector
////	}
//
//	return Vwind*1e-3;	//convert to km/s
//}
//
//
//
//ColVec Trajectory::CalcVinfECI(double timestep, ColVec State) {
//	ColVec Vinf(3,0.);
//	
//	double UTC_now = Initial_UTC + timestep/(24.*3600);
//	
//	//Unpack State vector (ECI)
//	ColVec r(3);			//in terms of [km]
//	r(0) = State(0);
//	r(1) = State(1);
//	r(2) = State(2);
//	ColVec v(3);			// [km/s]
//	v(0) = State(3);
//	v(1) = State(4);
//	v(2) = State(5);
//	
//	ColVec ECI__omega__ECEF(3,0.);
//	ECI__omega__ECEF(2) = rotEarthRad;
//	
//	switch (Wind_Option) {
//		case VACUUM: {
//			Vinf = v*(-1);  }
//			break;
//			
//		case SIMPLE_ATMOSPHERE: {
//			ColVec VatmECI(3,0.);
//			VatmECI = Cross_Vectors(ECI__omega__ECEF,r);	
//			Vinf = v*(-1) + VatmECI;	}
//			break;
//			
//			
//		case GRAM_WINDS: {
//			ColVec VwindUVW(3,0.);
//			
//			MatrixM ECI__C__ECEF(3,3,0.);
//			ECI__C__ECEF = ECEF_to_ECI_rotation(UTC_now);
//			
//			if (StateDot_Option == DEBRIS) {
//				gdlat = current_gdlat;
//				lon = current_lon;
//				VwindUVW = VwindCurveUVW(); }
//
//			else if (StateDot_Option == FIRST_STAGE) {
//				VwindUVW = VwindCurveUVW(r); }
//
//			//have to unpeel these multiplications otherwise boost complains
//			MatrixM UVW__C__ECEF(3,3,0.);
////			UVW__C__ECEF = prod(Rotation(1, -gdlat), prod(Rotation(1, PI/2) , prod(Rotation(3, PI/2) , Rotation(3, lon))));
//			UVW__C__ECEF = prod(Rotation(3, PI/2) , Rotation(3, lon));
//			UVW__C__ECEF = prod(Rotation(1, PI/2), UVW__C__ECEF);
//			UVW__C__ECEF = prod(Rotation(1, -gdlat), UVW__C__ECEF);
//			
////			cout << "gdlat = " << gdlat << endl;
////			exit(50);
//			
//			ColVec VwindECI(3,0.);
//			VwindECI = prod( trans (UVW__C__ECEF), VwindUVW );
//			VwindECI = prod(ECI__C__ECEF, VwindECI) + Cross_Vectors(ECI__omega__ECEF,r);
////			VwindECI = prod(ECI__C__ECEF, prod( trans (UVW__C__ECEF), VwindUVW ) ) + Cross_Vectors(ECI__omega__ECEF,r);
//			Vinf = v*(-1) + VwindECI;	}
//		break;
//
//			
//		default:
//			cout << "ERROR: CalcVinfECI never should have got here!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//			break;
//	}
//	
//	return Vinf;		// [km/s]
//}
//
//
//double Trajectory::DensityCurve(double altitude) {
//	double rho = -55;
//
//	
//	//rand_num->generate_number(sigma);
//	
//	switch (Density_Option) {
//		case VACUUM:
//			rho = 0.;
//			break;
//
//		case GRAM_DENSITY:
//			if ((altitude >= min_alt_rho) && (altitude <= max_alt_rho)) {
//				rho = gsl_spline_eval (RhoMean_spline, altitude, RhoMean_acc);
//			}
//			else if (altitude < min_alt_rho){
//				rho = gsl_spline_eval (RhoMean_spline, 0., RhoMean_acc); 
////				cout << "you may have crashed rho = " << rho << endl; 
//			}
//			else {			// way out in space
//				rho = 0; 
//			}
//			break;
//			
//		case CANTWELL_DENSITY:
//			if ((altitude >= 0)) {
//				rho = 1.225*exp(-755.30628334321*altitude/6371.0085);
//			}
//			else if (altitude < 0){
//				rho = 1.225*exp(-755.30628334321*0);
////				cout << "you may have crashed rho = " << rho << endl; 
//			}
//			else {			// way out in space
//				rho = 0; 
//			}			
//			break;
//			
//		case GRAM_DENSITY_UNCERTAIN:
//			if ((altitude >= min_alt_rho) && (altitude <= max_alt_rho)) {
//				
//				rho = gsl_spline_eval (RhoUncert_spline, altitude, RhoUncert_acc);
////				cout << "rho = " << rho << endl;
//				
//			}
//			else if (altitude < min_alt_rho){
//				rho = gsl_spline_eval (RhoUncert_spline, 0., RhoUncert_acc); 
//				//				cout << "you may have crashed rho = " << rho << endl; 
//			}
//			else {			// way out in space
//				rho = 0; 
//			}
//			break;
//			
//		default:
//			cout << "ERROR: You picked an invalid atmospheric density option!!! ~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//			break;
//	}
//	
//	return rho;		//units are [kg/m^3]
//}
//
//
//ColVec Trajectory::CalcAeroDrag(double timestep, ColVec State){
//	// BE EXTRA CAREFUL WITH UNITS IN HERE!!!
//
//	// Before entering this function, you should have already set the values for
//	//  * Aref
//	//  * current_alt
//	//  * current_gdlat
//	//  * current_lon
//	//  * current_radius_earth
//	//	* current_UTC
//
//	ColVec AeroDrag(3,0.);
//	
//	//Unpack State vector (ECI)
//	ColVec r(3);			//in terms of [km]
//	r(0) = State(0);
//	r(1) = State(1);
//	r(2) = State(2);
//	double mass = State(6);
//
////	double Cd = 0.2;
////	double Aref = 5;	//[m^2]
//	
//	//Aref: reference area of object.  For a rocket, this gets read in from the nominal trajectory file.
//		// For a piece of debris, this gets set in the debris initialization that happens for every piece.
//	
//	double earth_radius;
//	if (StateDot_Option == DEBRIS) {
////		double altitude = VecNorm(r) - earth_radius;
//		double rho = DensityCurve(current_alt);	//[kg/m^3]	 -- should return zero for vacuum
//
//		ColVec VinfECI(3,0.);
//		VinfECI = CalcVinfECI(timestep, State) * 1e3;	//converted to m/s
//
//		// Note: I didn't forget the |Vinf|^2, it's that one of them canceled out with Vinf/|Vinf|
//		AeroDrag =  VinfECI * (0.5 * rho * VecNorm(VinfECI) * Cd * Aref);	//[kg m/s^2]
//
//		//Convert into acceleration and m->km
//		AeroDrag = AeroDrag * (1e-3)*(1/mass);		// [km/s^2]
//
//	}
//	else if (StateDot_Option == FIRST_STAGE) {
//		//Throwing all of this inside an if-statement so that I can freely make changes for the debris case without breaking the first_stage one
//		//	Eventually will want this to be brought into agreement with debris case.
//		//  I think the only real difference is that don't need earth initial radius anymore.
//		earth_radius = initial_radius;
//
//		double altitude = VecNorm(r) - earth_radius;
//		double rho = DensityCurve(altitude);	//[kg/m^3]	 -- should return zero for vacuum
//
//		ColVec VinfECI(3,0.);
//		VinfECI = CalcVinfECI(timestep, State) * 1e3;	//converted to m/s
//
//		// Note: I didn't forget the |Vinf|^2, it's that one of them canceled out with Vinf/|Vinf|
//		AeroDrag =  VinfECI * (0.5 * rho * VecNorm(VinfECI) * Cd * Aref);	//[kg m/s^2]
//
//		//Convert into acceleration and m->km
//		AeroDrag = AeroDrag * (1e-3)*(1/mass);		// [km/s^2]
//	}
//    else if (StateDot_Option == WHOLE_ROCKET) {
//		//Throwing all of this inside an if-statement so that I can freely make changes for the debris case without breaking the first_stage one
//		//	Eventually will want this to be brought into agreement with debris case.
//		//  I think the only real difference is that don't need earth initial radius anymore.
//		earth_radius = initial_radius;
//        
//		double altitude = VecNorm(r) - earth_radius;
//		double rho = DensityCurve(altitude);	//[kg/m^3]	 -- should return zero for vacuum
//        
//		ColVec VinfECI(3,0.);
//		VinfECI = CalcVinfECI(timestep, State) * 1e3;	//converted to m/s
//        
//		// Note: I didn't forget the |Vinf|^2, it's that one of them canceled out with Vinf/|Vinf|
//		AeroDrag =  VinfECI * (0.5 * rho * VecNorm(VinfECI) * Cd * Aref);	//[kg m/s^2]
//        
//		//Convert into acceleration and m->km
//		AeroDrag = AeroDrag * (1e-3)*(1/mass);		// [km/s^2]
//	}
//    else if (StateDot_Option == SUBORBITAL) {
//		//Throwing all of this inside an if-statement so that I can freely make changes for the debris case without breaking the first_stage one
//		//	Eventually will want this to be brought into agreement with debris case.
//		//  I think the only real difference is that don't need earth initial radius anymore.
//		earth_radius = initial_radius;
//        
//		double altitude = VecNorm(r) - earth_radius;
//		double rho = DensityCurve(altitude);	//[kg/m^3]	 -- should return zero for vacuum
//        
//		ColVec VinfECI(3,0.);
//		VinfECI = CalcVinfECI(timestep, State) * 1e3;	//converted to m/s
//        
//		// Note: I didn't forget the |Vinf|^2, it's that one of them canceled out with Vinf/|Vinf|
//		AeroDrag =  VinfECI * (0.5 * rho * VecNorm(VinfECI) * Cd * Aref);	//[kg m/s^2]
//        
//		//Convert into acceleration and m->km
//		AeroDrag = AeroDrag * (1e-3)*(1/mass);		// [km/s^2]
//        
////        cout << "t = " << timestep << "  alt = " << altitude << "   vel = " << VecNorm(VinfECI) << endl;
//	}
//	else {
//		cout << "CalcAeroDrag should never have got here~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//
//	
//	return AeroDrag;
//}














//
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Possible States to Integrate and Integration Methods
//
//ColVec Trajectory::Statedot_first_stage(double timestep, ColVec State) {
//	
//	//Unpack State Vector
//	ColVec r(3);			//in terms of [km]
//	r(0) = State(0);
//	r(1) = State(1);
//	r(2) = State(2);
//	ColVec v(3);			// [km/s]
//	v(0) = State(3);
//	v(1) = State(4);
//	v(2) = State(5);
//	double mass = State(6);	//[kg]
//	
//	// Acceleration due to center-pointing gravity
//	ColVec acc_grav(3);
//	acc_grav = r * (-mu_earth / pow(VecNorm(r),3));
//	   
//    
//	// Force due to thrust
//	ColVec ThrustForce(3,0.);
//	ThrustForce = GetThrustForce(timestep) + DeltaForce;
//	double m_dot = -VecNorm(ThrustForce)/(Isp*9.8);	//[kg/s]
//    
//
//    
//    
//    
//    
//	
//	// From force, can calc acceleration due to thrust
//	ColVec acc_thrust(3,0.);
//	acc_thrust = ThrustForce*(1/(mass))*(1e-3);	//1e3 is to convert thrust into km/s^2
//	
//	// Force due to aerodynamic drag	-- Should return zero for vacuum
//	ColVec acc_drag(3,0.);
//	acc_drag = CalcAeroDrag(timestep, State);
//	
//	// Total acceleration is sum of previous acceleration vectors
//	ColVec acc_total(3);
//	acc_total = acc_grav + acc_thrust + acc_drag;
//	
//	//	cout << "grav = " << acc_grav;
//	//	cout << "thrust = " << acc_thrust;
//	//	cout << "drag = " << acc_drag;
//	
//	// Pack up
//	ColVec Statedot(state_vec_size,0.0);
//	
//	Statedot(0) = State(3);
//	Statedot(1) = State(4);
//	Statedot(2) = State(5);
//	
//	Statedot(3) = acc_total(0);
//	Statedot(4) = acc_total(1);
//	Statedot(5) = acc_total(2);
//	
//	Statedot(6) = m_dot;
//	
//	return Statedot;
//}
//
//ColVec Trajectory::Statedot_debris(double timestep, ColVec State) {
//	
//	//Unpack State Vector
//	ColVec r(3);			//in terms of [km]
//	r(0) = State(0);
//	r(1) = State(1);
//	r(2) = State(2);
//	ColVec v(3);			// [km/s]
//	v(0) = State(3);
//	v(1) = State(4);
//	v(2) = State(5);
////	double mass = State(6);	//[kg]
//	
//	//Calc some useful properties of the given state
//	current_UTC = Initial_UTC + timestep / (3600*24);
//	current_ECEF = prod( trans(ECEF_to_ECI_rotation(current_UTC)), r );
////	cout << "current_ECEF = " << current_ECEF << endl;
//	ColVec Geod = ECEF_To_Geodetic(current_ECEF);
////	cout << "Gedo = " << Geod << endl;
//	current_gdlat = Geod(0);
//	current_lon = Geod(1);
//	current_alt = Geod(2);
//	current_radius_earth = get_local_earth_radius(gdlat);
//
//
//	// Acceleration due to center-pointing gravity
//	ColVec acc_grav(3);
//	acc_grav = r * (-mu_earth / pow(VecNorm(r),3));
//	
//	// Force due to thrust
//	double m_dot = 0;	//[kg/s]
//	
//	// Force due to aerodynamic drag	-- Should return zero for vacuum
//	ColVec acc_drag(3,0.);
//	acc_drag = CalcAeroDrag(timestep, State);
//	
//	// Total acceleration is sum of previous acceleration vectors
//	ColVec acc_total(3);
//	acc_total = acc_grav + acc_drag;
//	
////	cout << "mag of acc_total = " << VecNorm(acc_total) << endl;
////	cout << "mag of r_dot = " << VecNorm(v) << endl;
//	
//	// Pack up
//	ColVec Statedot(state_vec_size,0.0);
//	
//	Statedot(0) = State(3);
//	Statedot(1) = State(4);
//	Statedot(2) = State(5);
//	
//	Statedot(3) = acc_total(0);
//	Statedot(4) = acc_total(1);
//	Statedot(5) = acc_total(2);
//	
//	Statedot(6) = m_dot;
//	
////	cout << "Statedot = " << Statedot << endl;
//	
//	return Statedot;
//}
//
//ColVec Trajectory::Statedot_suborbital(double timestep, ColVec State) {
//	
//	//Unpack State Vector
//	ColVec r(3);			//in terms of [km]
//	r(0) = State(0);
//	r(1) = State(1);
//	r(2) = State(2);
//	ColVec v(3);			// [km/s]
//	v(0) = State(3);
//	v(1) = State(4);
//	v(2) = State(5);
//	double mass = State(6);	//[kg]
//	
//	// Acceleration due to center-pointing gravity
//	ColVec acc_grav(3);
//	acc_grav = r * (-mu_earth / pow(VecNorm(r),3));
//    
//	
//	// Force due to thrust
////	ColVec ThrustForce(3,0.);
//    // TODO: HERE
////	ThrustForce = GetThrustForce(timestep);
//    
//    ColVec acc_body(3,0.);
//    acc_body = GetBodyAcc(timestep);
//    double body_acc_norm = VecNorm(acc_body);   //Save the total body acceleration for mdot
////    acc_body(2) = 0.;   //Per my convention, throw out the thrust in the z-dir
//    
//    ColVec bodyAngles(3,0.);    // (pitchAngle, omega, omegaDot)
//    bodyAngles = GetBodyAngles(timestep);
//    
////    cout << "DEBUGGIN pitchAngle!!!!~~~~~~~~~~~~~~~~\n";
//    //BUG: Seems like pitchAngle = 0 shoots straight up, 90 runs it due east (into the ground), -90 goes west
//    double pitchAngle = bodyAngles(0);
////    double omega = bodyAngles(1);
////    double omegaDot = bodyAngles(2);
//    
//    ColVec eci_V_s__eci(3,0.);
//    eci_V_s__eci = v;
//    
//    ColVec o_R_s__eci(3,0.);
//    o_R_s__eci = r;
//    
//    ColVec eci_omega_ecef__eci(3);
//	eci_omega_ecef__eci(0) = 0;
//	eci_omega_ecef__eci(1) = 0;
//	eci_omega_ecef__eci(2) = rotEarthRad;
//    
//    //Rotation matrix for converting from ecef to eci
//	double UTC = Initial_UTC;
//	MatrixM eci__C__ecef(3,3);
//	eci__C__ecef = ECEF_to_ECI_rotation(UTC);
//    
//    MatrixM ecef__C__eci(3,3);
//    ecef__C__eci = trans(eci__C__ecef);
//    
//    ColVec o_R_s__ecef(3,0.);
//    o_R_s__ecef = prod(ecef__C__eci, o_R_s__eci);
//    
//    ColVec LatLonAlt(3,0.);
//    LatLonAlt = ECEF_To_Geodetic(o_R_s__ecef);
//    
//    //Rotate to ECEF (which is fixed wrt to local)
//    MatrixM local__C__ecef (3,3,0.);
//    local__C__ecef = prod(Rotation(2, -LatLonAlt(0)), Rotation(3, LatLonAlt(1)));
//    
//    //This is good
//    ColVec ecef_V_s__eci(3,0.);
//    ecef_V_s__eci = eci_V_s__eci - Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci); //CrossVec is good
//    
//    
////    cout << "DEBUG\n eci_V_s__eci = " << eci_V_s__eci << " ,  CrossVec = " << Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci) << endl;
//    
//    ColVec ecef_V_s__local(3,0.);
//    ecef_V_s__local = prod(ecef__C__eci, ecef_V_s__eci);    //This step is confirmed to be ecef_V_s__ecef
//    ecef_V_s__local = prod(local__C__ecef, ecef_V_s__local);
//    
//    double headingAngle = atan2(ecef_V_s__local(2), ecef_V_s__local(1));
////    double headingAngle = 90.*PI/180;   //Currently doesn't do ANYTHING?!
//    
//    //Now rotate body frame by -headingAngle
//    MatrixM body__C__local (3,3,0.);
////    body__C__local = prod(Rotation(3, -pitchAngle * PI/180), Rotation(1, headingAngle));
//    body__C__local = prod(Rotation(2, -pitchAngle * PI/180), Rotation(3, headingAngle));
//    
//    ColVec local_A_s__localTEMP(3,0.);
//    local_A_s__localTEMP = prod(trans(body__C__local), acc_body);
//    
//    ColVec local_A_s__local(3,0.);
//    local_A_s__local(1) = local_A_s__localTEMP(0);
//    local_A_s__local(2) = local_A_s__localTEMP(1);
//    local_A_s__local(0) = local_A_s__localTEMP(2);
//    
////    cout << "local_A_s__local = " << local_A_s__local << endl;
//    
////    cout << "latlonalt = " << LatLonAlt << endl;
//    
//    ColVec local_A_s__eci(3,0.);
//    local_A_s__eci = prod( trans(local__C__ecef), local_A_s__local);  //local_A_s__eci is currently set to local_A_s__ecef
//    local_A_s__eci = prod( eci__C__ecef, local_A_s__eci);               //now local_A_s__eci is correct
//    
//    ColVec coriolis__eci(3,0.);
//    coriolis__eci = 2*Cross_Vectors(eci_omega_ecef__eci, ecef_V_s__eci);
//    
//    ColVec centrifugal__eci(3,0.);
//    centrifugal__eci = Cross_Vectors(eci_omega_ecef__eci, Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci));  //might throw an error
//    
//    ColVec eci_A_s__eci(3,0.);
//    eci_A_s__eci = (local_A_s__eci*(1e-3)*9.8 + coriolis__eci + centrifugal__eci); //(9.8 to convert g's to m/s^2) (1e3 to convert acceleration into km/s^2)
//    
////    cout << "eci_A_s__eci = " << VecNorm(local_A_s__eci)*(1e-3)*9.8 << " + " << VecNorm(coriolis__eci)
////    << " + " <<  VecNorm(centrifugal__eci) << " , headingAngle = "  << headingAngle *180/PI << " , pitchAngle = " << pitchAngle << endl;
//
//    
////    eci_V_s__eci = prod(eci__C__ecef,ecef_V_s__ecef) + Cross_Vectors(eci_omega_ecef__eci, o_R_s__eci);
//
//    
//    
//    
//    
//    
//    
//    
//    
//    
////	double m_dot = -VecNorm(ThrustForce)/(Isp*9.8);	//[kg/s]
//	double m_dot = -body_acc_norm*mass/(Isp);	//[kg/s] note that acc_body is measured in g's
//	
////	// From force, can calc acceleration due to thrust
////	ColVec acc_thrust(3,0.);
////	acc_thrust = acc_body*(1e-3)*9.8;	//(9.8 to convert g's to m/s^2) (1e3 to convert acceleration into km/s^2)
//	
//	// Force due to aerodynamic drag	-- Should return zero for vacuum
//	ColVec acc_drag(3,0.);
//	acc_drag = CalcAeroDrag(timestep, State);
//	
//	// Total acceleration is sum of previous acceleration vectors
//	ColVec acc_total(3);
//	acc_total = acc_grav + eci_A_s__eci + acc_drag;
//	
////    cout << "total = " << VecNorm(acc_grav) << " + " << VecNorm(eci_A_s__eci) << " + " <<  VecNorm(acc_drag) << " = " << VecNorm(acc_total) << endl;
//	
//	// Pack up
//	ColVec Statedot(state_vec_size,0.0);
//	
//	Statedot(0) = State(3);
//	Statedot(1) = State(4);
//	Statedot(2) = State(5);
//	
//	Statedot(3) = acc_total(0);
//	Statedot(4) = acc_total(1);
//	Statedot(5) = acc_total(2);
//	
//	Statedot(6) = m_dot;
//	
//	return Statedot;
//}
//
//
//
//
////Note: the incoming params pointer is actually the 'this' pointer.  That's because we're trying to call some
//// class member functions from within a static function, so we need an object in the class from which to call them.
//int Trajectory::func (double t, const double y[], double f[], void *params) {	
//	Trajectory *current_obj = (Trajectory *)params;
//	
//	//load up the y vector into a colvec
//	ColVec current_state(state_vec_size);
//	for (int i = 0; i < state_vec_size; i++) {
//		current_state(i) = y[i]; }
//
//	
//	ColVec statedot(state_vec_size,0.);	
//	switch (current_obj->StateDot_Option) {
//        case WHOLE_ROCKET:
//			statedot = current_obj->Statedot_first_stage(t,current_state);
//			break;
//            
//		case FIRST_STAGE:
//			statedot = current_obj->Statedot_first_stage(t,current_state);
//			break;
//			
//		case DEBRIS:
//			statedot = current_obj->Statedot_debris(t,current_state);
//			break;
//            
//        case SUBORBITAL:
//            statedot = current_obj->Statedot_suborbital(t,current_state);
//            break;
//
//		default:
//			cout << "ERROR: INVALID STATEDOT!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
//			break;
//	}
//	
//	
//	// unpack the ColVec
//	f[0] = statedot(0);
//	f[1] = statedot(1);
//	f[2] = statedot(2);
//	f[3] = statedot(3);
//	f[4] = statedot(4);
//	f[5] = statedot(5);
//	f[6] = statedot(6);
//		
//	return GSL_SUCCESS;
//}
//
//
//
//
//int Trajectory::gsl_integrate_stages (ColVec State0, int run_number) {
// 
//    // Set up the gsl integration
//    gsl_odeiv2_system sys = {Trajectory::func, NULL, state_vec_size, this};
//    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
//    
//    // Keep track of where the integration ends
//    int finalStateIndex = ((INTxx) RocketStateTimes.size())-1;
//
//    // Different integrators
//    // gsl_odeiv2_step_rk2
//    // *gsl_odeiv2_step_rkf45
//    
//    double y[] = {State0(0), State0(1), State0(2), State0(3), State0(4), State0(5), 0};
//    int tcount = 1;  //skip the t=0 point //this counts the state times
//
//    for (int stageNum = 0; stageNum < NumStages; stageNum++) {
//        
//        // This is the time that the integration starts at (gets updated for each timestep)
//        double t0 = StageTimes[stageNum];
//        
//        // The initial weight for the rocket at the beginning of each stage
//        double StageMass = PayloadMass;
//        for (int jx = stageNum; jx < NumStages; jx++) {
//            StageMass += StagePropMass[jx] + StageStructMass[jx]; }
//        
//        //initialize y_vec
//        y[6] = StageMass;
//        
////        double y[] = {State0(0), State0(1), State0(2), State0(3), State0(4), State0(5), StageMass};
//        std::vector <double> yvecstore (y, y + state_vec_size); //Copies y into yvecstore (need to add array size to initialize from an array)
//        if ( stageNum == 0) {
//            State_Vector_Storage_Vec[run_number][0] = yvecstore; }
//        
//        // Though poorly named, this variable specifies the index up to which the propagation has been computed
////        debris_state_index = StateTimesLength-1;
//        
//        cout << "numtimestep = " << RocketNumTimeSteps[stageNum] << endl;
//        
//        for (int ixx = 1; ixx < RocketNumTimeSteps[stageNum]; ixx++)
//        {
////            cout << "ixx = " << ixx << endl;
////            static int tcount = 1;  //skip the t=0 point //this counts the state times
////            cout << "tcount = " << tcount << endl;
//
//            double ti = RocketStateTimes[tcount];
//            //		cout << "i = " << i << endl;
//            int status = gsl_odeiv2_driver_apply (d, &t0, ti, y);	//BEWARE!!!  Function changes t0
//            
//            double current_alt = pow(y[0]*y[0] + y[1]*y[1] + y[2]*y[2],0.5) - initial_radius;	//[km]
//            double altitude_limit = 2000; //18.288;	//[km]
//            
//            //store the current solution
//            std::vector <double> yvecstore (y, y + state_vec_size);
//            
////            cout << "Current Mass = " << y[6] << endl;
//            
//            State_Vector_Storage_Vec[run_number][tcount] = yvecstore;
//            
////            if (State_Vector_Storage_Vec[0][0][6] != 3960.) {
////                cout << "enter debug" << endl;
////            }
//            
//            if (current_alt > altitude_limit){
//                cout << "rocket has left upper limit of " << altitude_limit << "km, at time = " << ti << " seconds" << endl;
////                debris_state_index = tcount;
//                finalStateIndex = tcount;
//                break; }
//            else if (current_alt < -1) {
//                cerr << "You have probably crashed!  seconds after launch = " << ti << endl;
////                debris_state_index = tcount;
//                finalStateIndex = tcount;
//                break; }
//            else if ((StageMass - y[6]) > StagePropMass[stageNum]) {
//                cerr << "You ran out of propellant!!!!   seconds after launch = " << ti << endl;
//                finalStateIndex = tcount;
////                break;
//            }
//            else if (status != GSL_SUCCESS)
//            {
//                printf ("error, return value=%s\n", gsl_strerror(status));
//                break; }
//            
//            finalStateIndex = tcount;
//
//            tcount++;
//        }
//        
//        
//        // Start propagating the staged masses down ----------
//        // Initialize the state vector storage
////        FallingStageStateVectorStorage[run_number][0] =
//        
//        // If statement that turns this on when we want first stage only!!!
//        break;
//        
//    }
//    
//    //Resize the vector if needed
//    State_Vector_Storage_Vec[run_number].resize(finalStateIndex+1);
//    
//    
////    if (stageNum == 0) {
////        y[6] = StageStructMass[0];
////        FallingStageStateVectorStorage[run_number] = Propagate_Stage_Down_To_Ground(y, StageTimes[1], 0);
////        cout << "exited the new func"   << endl;
////    }
//    
//
//	
//	//UNCOMMENT THIS WHEN DONE TRYING FIXED TIMESTEP!!!
//	// Code seems to be a little bit faster with this commented out, if it doesn't impact performance
//	//   then move this command someplace else.
//	gsl_odeiv2_driver_free (d);
//    
//	return 0;
//}
//
//
//
//
//
//
//
//int Trajectory::gsl_integrate_first_stage (ColVec State0, int run_number) {
//		
//	double t0 = time_start;
//		
//	// Different integrators
//	// gsl_odeiv2_step_rk2
//	// *gsl_odeiv2_step_rkf45
//	
//	gsl_odeiv2_system sys = {Trajectory::func, NULL, state_vec_size, this};
//	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
//	
//	//initialize y_vec
//	double y[] = {State0(0), State0(1), State0(2), State0(3), State0(4), State0(5), State0(6)};
//	std::vector <double> yvecstore (y, y + state_vec_size);
//	State_Vector_Storage_Vec[run_number][0] = yvecstore;
//
////	int debris_state_index = rocket_num_time_steps-1;
//	
//    int rocket_num_time_steps = (INTxx) State_Vector_Storage_Vec[run_number].size();
//	cout << "numtimestep = " << rocket_num_time_steps << endl;
//
//	for (int i = 1; i < rocket_num_time_steps; i++)
//	{
//		double ti = RocketStateTimes[i];
////		cout << "i = " << i << endl;
//		int status = gsl_odeiv2_driver_apply (d, &t0, ti, y);	//BEWARE!!!  Function changes t0
//		
//		double current_alt = pow(y[0]*y[0] + y[1]*y[1] + y[2]*y[2],0.5) - initial_radius;	//[km]
//		double altitude_limit = 200; //18.288;	//[km]
//
//		//store the current solution 
//		std::vector <double> yvecstore (y, y + state_vec_size);
//		
//		State_Vector_Storage_Vec[run_number][i] = yvecstore;
//		
//		if (current_alt > altitude_limit){
//			cout << "rocket has left upper limit at time = " << ti << " seconds" << endl;
//			debris_state_index = i;
//			break; }
//		else if (current_alt < -1) {
//			cerr << "You have probably crashed!  seconds after launch = " << ti << endl;
//			debris_state_index = i;
//			break; }
//		else if (status != GSL_SUCCESS)
//     	{
//			printf ("error, return value=%s\n", gsl_strerror(status));
//			break;
//     	}
//
//
//	}
//		
////	cout << "debris_state = {";
////	for (int i = 0; i < state_vec_size; i++) {
////		cout << std::setprecision(16) << State_Vector_Storage_Vec[run_number][debris_state_index][i] << ", "; }
////	cout << Initial_UTC + RocketStateTimes[debris_state_index]/(24*3600) << "};" << endl;
//	
//	//UNCOMMENT THIS WHEN DONE TRYING FIXED TIMESTEP!!!
//	// Code seems to be a little bit faster with this commented out, if it doesn't impact performance
//	//   then move this command someplace else.
//	gsl_odeiv2_driver_free (d);
//
//	return 0;
//}
//
//// I think this function is nearly identical to the stages one
//int Trajectory::gsl_integrate_suborbital (ColVec State0, int run_number) {
//    // TODO: This function is exactly the same as first_stage.  Clean this up via merge and deletion
//	double t0 = time_start;
//    
//    int StateTimesLength = (INTxx) State_Vector_Storage_Vec[run_number].size();
//	cout << "numtimestep = " << StateTimesLength << endl;
//    
//    // Keep track of where the integration ends
//    int finalStateIndex = StateTimesLength-1;
//    
//	// Different integrators
//	// gsl_odeiv2_step_rk2
//	// *gsl_odeiv2_step_rkf45
//	
//	gsl_odeiv2_system sys = {Trajectory::func, NULL, state_vec_size, this};
//	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
//	
//	//initialize y_vec
//	double y[] = {State0(0), State0(1), State0(2), State0(3), State0(4), State0(5), State0(6)};
//	std::vector <double> yvecstore (y, y + state_vec_size);
//	State_Vector_Storage_Vec[run_number][0] = yvecstore;
//    
////	debris_state_index = rocket_num_time_steps-1;
//    int IX = 0;
//    for (int i = 1; i < StateTimesLength; i++)
//    {
//        IX++;
//        
//        double ti = SuborbitalStateTimes[i];
//        cout << "i = " << i << endl;
//        int status = gsl_odeiv2_driver_apply (d, &t0, ti, y);	//BEWARE!!!  Function changes t0
//        
//        double current_alt = pow(y[0]*y[0] + y[1]*y[1] + y[2]*y[2],0.5) - initial_radius;	//[km]
//        double altitude_limit = 200; //18.288;	//[km]
//        
//        //store the current solution
//        std::vector <double> yvecstore (y, y + state_vec_size);
//        
//        State_Vector_Storage_Vec[run_number][i] = yvecstore;
//        
//        if (current_alt > altitude_limit){
//            cout << "rocket has left upper limit at time = " << ti << " seconds" << endl;
//            debris_state_index = i;
//            finalStateIndex = i;
//            break; }
//        else if (current_alt < -1) {
//            cerr << "You have probably crashed!  seconds after launch = " << ti << endl;
//            debris_state_index = i;
//            finalStateIndex = i;
//            break; }
//        else if (status != GSL_SUCCESS)
//        {
//            printf ("error, return value=%s\n", gsl_strerror(status));
//            break;
//        }
//        
//        
//    }
//    
////    cout << "debris_state = {";
////    for (int i = 0; i < state_vec_size; i++) {
////        cout << std::setprecision(16) << State_Vector_Storage_Vec[run_number][debris_state_index][i] << ", "; }
////    cout << Initial_UTC + SuborbitalStateTimes[debris_state_index]/(24*3600) << "};" << endl;
//    
//    //Resize the vector if needed
////    State_Vector_Storage_Vec[run_number].resize(finalStateIndex+1);
//    State_Vector_Storage_Vec[run_number].resize(IX);
//	
//	//UNCOMMENT THIS WHEN DONE TRYING FIXED TIMESTEP!!!
//	// Code seems to be a little bit faster with this commented out, if it doesn't impact performance
//	//   then move this command someplace else.
//	gsl_odeiv2_driver_free (d);
//    
//	return IX;
//}
//
//
//
//
//int Trajectory::gsl_integrate_debris (ColVec State0, double *StateTimes) {
//	
//	double t0 = time_start;     //Is there any time that I'll ever need a time_start other than zero???
//	
//	gsl_odeiv2_system sys = {Trajectory::func, NULL, state_vec_size, this};
//	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
//	
//	//initialize y_vec
//	double y[] = {State0(0), State0(1), State0(2), State0(3), State0(4), State0(5), State0(6)};
//	std::vector <double> yvecstore (y, y + state_vec_size);
//	State_Vector_Buffer_Vec[0] = yvecstore;
//	
//	// Going to take the gdlat of the initial state vector to find the approximate local earth radius
//	// Think about this a little more when you have the time
//	ColVec GeoCoords = ECEF_To_Geodetic(State0);	//State0 is technically in ECI, but that won't affect alt calculation
////	double earth_radius = 6371.0008;
////	double earth_radius = get_local_earth_radius(GeoCoords(0));
//	
////	double current_alt = pow(y[0]*y[0] + y[1]*y[1] + y[2]*y[2],0.5) - earth_radius;	//[km]
//	double current_alt = GeoCoords(2);	//[km]
//	
////	cout << "current time = " << t0 << endl;
////	cout << "current_alt = " << current_alt << endl;
//	
//	int IX = 0;
//	while (current_alt > 0)
//	{
//		IX++;						// So we're actually starting at i=1
//		
//		double ti = StateTimes[IX];
//		
//		int status = gsl_odeiv2_driver_apply (d, &t0, ti, y);	//BEWARE!!!  Function changes t0
//		
//		GeoCoords = ECEF_To_Geodetic(y[0],y[1],y[2]);
//		current_alt = GeoCoords(2);
//        
////		current_alt = pow(y[0]*y[0] + y[1]*y[1] + y[2]*y[2],0.5) - earth_radius;	//[km]
//
//		if (status != GSL_SUCCESS)
//     	{
//			printf ("error, return value=%s\n", gsl_strerror(status));
//			break;
//     	}
//		
//		//store the current solution 
//		std::vector <double> yvecstore (y, y + state_vec_size);
//		State_Vector_Buffer_Vec[IX] = yvecstore;
//	}
//	
//    
////	cout << "Your debris hit the earth at time = " << StateTimes[IX] << endl;
////	cout << "current_alt = " << current_alt << "   at timestep = " << IX << endl;
//
//	
//	// Loop did not continue, that means IX entry is the last one
//	//  When converting to altitude, should make this (and probable subsequent entries) be zero
//	
//	
//	//UNCOMMENT THIS WHEN DONE TRYING FIXED TIMESTEP!!!
//	// Code seems to be a little bit faster with this commented out, if it doesn't impact performance
//	//   then move this command someplace else.
//	gsl_odeiv2_driver_free (d);
//	
//	return IX;
//}
//















//
//
//void Trajectory::write_LRHC_points_file(char *outFileName, double tstepMin) {
//
////    write_points_file(outFileName, tstepMin, WHOLE_ROCKET);
//    write_points_file(outFileName, tstepMin, StateDot_Option);
//    return;
//}
//
//void Trajectory::write_stage_down_points_file(char *outFileName, double tstepMin) {
//    
//    write_points_file(outFileName, tstepMin, STAGE_DOWN);
//    return;
//}
//
//void Trajectory::write_points_file(char *outFileName, double tstepMin, int StateDot_In) {
//    // This function writes the LRHC ellipse around the nominal trajectory to a points file so that it may be wrapped up into a footprint
//    // Here, we bin the points according to time
//    //  ASSUMES: Traditional rocket that goes up and never re-enters NAS; will break if you use this for something that re-enters
//    //  ASSUMES: Your first timestep has the rocket within the NAS
//	
//	// Do some error checking
//	ofstream outfile;
//	outfile.open(outFileName,ios::out | ios::binary);
//	
//	if (outfile.bad()) {
//		cerr << "STAGE DOWN FILE FAILED TO OPEN FOR OUTPUT!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//	
//	// Things we'll need
//	Point temp_pt;
//	std::vector<std::vector<Point> > total_points_at;	//[timestep][points]
//    
//    // We'll want to output this, so save it
//    double startUTC = Initial_UTC;
//    
//    int numRuns = num_per_batch;
//
//    std::vector<std::vector<std::vector<double> > > LatLonAltStorage;
//    
//    
//    double timeInFlight;
//    int old_StateDot_Option = StateDot_Option;
//    
//    // Figure out which statetimes we need
//    vector <double> StateTimes;
//    if (StateDot_In == FIRST_STAGE){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_In == WHOLE_ROCKET){
//        StateTimes = RocketStateTimes;
//        
//        // This function is going to only consider times up to the first staging event
//        timeInFlight = (StageTimes[1] - StateTimes[0])/60.;	//minutes
//        
//        StateDot_Option = WHOLE_ROCKET; //Change to desired option
//        LatLonAltStorage = TransformToLatLonAlt(numRuns);
//
//    
//    }
//    else if (StateDot_In == SUBORBITAL) {
//        
//        StateTimes = SuborbitalStateTimes;
//        timeInFlight = (StateTimes.back() - StateTimes[0])/60.;	//minutes
//        LatLonAltStorage = TransformToLatLonAlt(numRuns);
//        
//        cout << "StateTimes.back() = " << StateTimes.back() << endl;
//        cout << "StateTimes[0] = " << StateTimes[0] << endl;
//        cout << "timeInFlight = " << timeInFlight << endl;
//
//        cout << "this function doesn't handle suborbital case yet, exiting 87234" << endl;
////        exit(23);
//      }
//    else if (StateDot_In == STAGE_DOWN) {
//        StateTimes = StageDownStateTimes;
//        timeInFlight = (StateTimes.back() - StateTimes[0])/60.;	//minutes
//        cout << "StateTimes.back() = " << StateTimes.back() << endl;
//        cout << "StateTimes[0] = " << StateTimes[0] << endl;
//        cout << "timeInFlight = " << timeInFlight << endl;
//        
//        StateDot_Option = STAGE_DOWN; //Change to desired option
//        LatLonAltStorage = TransformToLatLonAlt(numRuns);
//
//
//        
//    }
//    StateDot_Option = old_StateDot_Option;  //Change back
//
//    
//    
//
//    
//    //Allocate the vector of points that we're going to write to file
//    
//    int time_steps_out = (INTxx) ceil(timeInFlight/tstepMin);
//    total_points_at.assign(time_steps_out,std::vector<Point>());
//    int num_range = time_steps_out;
//    
//
//    
//    Point tempPoint;
//    
//    double altitudeLimit = 18.288;      // [km] This is the top of the NAS
//    
//    for (int runNumber = 0; runNumber < numRuns; runNumber++){
//        int numPtsHere = (INTxx) LatLonAltStorage[runNumber].size();
//        for (int t = 0; t < numPtsHere; t++){
//            double curAlt = LatLonAltStorage[runNumber][t][2];
//            if (curAlt <= altitudeLimit) {
//                int curTimeIX = (INTxx) floor((StateTimes[t] - StateTimes[0])/(60*tstepMin));
//                
//                // Kludge.  It's possible that timeInFlight/tstepMins is an integer, thus ceil doesn't round up, thus you can hit the upper limit exactly
//                //      and thus try accessing memory that's out of range.  This if statement catches that case, though it's an ugly solution.
//                if (curTimeIX < time_steps_out) {
//                    double gdlat = LatLonAltStorage[runNumber][t][0];
//                    double lon = LatLonAltStorage[runNumber][t][1];
//                    double rLocal = get_local_earth_radius(gdlat);
//                    
//                    // Anticipate some problems here
//                    if (curAlt == 0) {
//                        curAlt = 0.001; //Set at 1m above ground to avoid division by zero
//                    }
//                    else if (curAlt < 0) {
//                        cout << "ERROR: should not have al0 negative in here, come fix this!!!!~~~~~~~~~" << endl;
//                    }
//                    
////                    // Location of nominal trajectory on the ground
////                    double curX = R_equator * LatLonAltStorage[runNumber][t][1];  // R_equator * longitude_angle
////                    double curY = rLocal * gdlat;                // Local Radius of Earth * gdlat
//                    
//                    tempPoint.set_Point(gdlat, lon, curAlt, rLocal);
//                    
////                    tempPoint.set_x(curX);
////                    tempPoint.set_y(curY);
//                    
////                    tempPoint.set_z(curAlt);
////                    tempPoint.set_R_local(rLocal);
////                    tempPoint.set_gdlat(gdlat);
////                    tempPoint.set_lon(LatLonAltStorage[runNumber][t][1]);
//                    
//                    total_points_at[curTimeIX].push_back(tempPoint);
//                }
//            }
//        } }
//    
//	
//	// ------------ Begin writing to file
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
//	
//	outfile.close();
//	
//	// Used for debugging purposes only
//	bool debug = false;
//	if (debug) {
//		string debugFileName("GeneratedFiles/LRHC_PLAINTEXT_DEBUG.txt");
//		ofstream debugfile;
//		debugfile.open(debugFileName.c_str(), ios::out);
//		
//		debugfile << num_range << endl << startUTC << endl << tstepMin << endl;
//		
//		int num_points_here;
//		for (int t = 0; t < num_range; t++) {
//			num_points_here = (INTxx) total_points_at[t].size();
//			for (int ix = 0; ix < num_points_here; ix++) {
//				debugfile << total_points_at[t][ix] << endl; }
//			debugfile << "\n\n\n\n\n\n\n\n\n\n\n\n\n"; }
//		
//		debugfile.close();
//	}
//	// End of debugging code
//	
//	
//	return;
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
////double Trajectory::AtWhatTimeCross(double thresholdKm){
////    
////}
//
//
//
//
//double Trajectory::HowLongUntilEnterNAS() {
//    int oldStateDot = StateDot_Option;
//    StateDot_Option = DEBRIS;
//    
//    int numPieces = DebrisLatLonAltStorage.size();
//
//    
//    std::vector<std::vector<std::vector<double> > > LatLonAltStorage;
//    LatLonAltStorage = TransformToLatLonAlt(numPieces);
//    
//    int shortestTime = 10000;   //timesteps, not seconds
//    double NAS = 18.288;    //height of NAS in km
//    
//    for (int px = 0; px < numPieces; px++) {
//        int numTSteps = LatLonAltStorage[px].size();
//
//        for (int ix = 0; ix < numTSteps; ix++){
//            if ((LatLonAltStorage[px][ix][2] < NAS) && (shortestTime > ix)){
//                shortestTime = ix;
//                cout << "shortestTime = " << shortestTime << endl;
//
//            }
////            cout << "Alt = " << LatLonAltStorage[0][ix][2] << endl;
//
//            
//        }
//    }
//    
//    cout << "numPieces = " << numPieces << endl;
//    cout << "shortest time in minutes = " << DebrisStateTimes[shortestTime]/60. << endl;
//    
//    
//    StateDot_Option = oldStateDot;
//    return 6.;
//    
//} 
//
//
//
///*! MonteCarloDebris: Currently, this function
// *  * Inputs a pre-existing footprint object and updates it
// *  *  * Might be a nice idea to have option allowing the footprint object to be updated
// *  *  * Would also be nice to pass in the state-vectors you want to integrate, rather than figuring them out in this function
// *  * propagates debris to the ground starting from altitude_we_explode_until_km at every tstep = tstep - backwardsStep
// *  * at each considered timestep (of initial explosion) APPENDS the generated points into the existing envelope
// *  * when finished, stores the footprint as a vector in a binary file via store_footprint_as_vector
// */
//void Trajectory::MonteCarloDebris(Footprint3D &my_footprint, double time_explode_until_in){
//    // By this point, the parameters of the architecture should be already defined.  Needed params are:
//    //    double binSizeKm = 5.;  //Actually, everything in Trajectory should be agnostic to this one
//    
//    double tstepMinutes = my_footprint.getDeltaT();
//    double launchLat = my_footprint.getLaunchLat();
//    double launchLon = my_footprint.getLaunchLon();
//    double launchAzimuth = my_footprint.getLaunchAzimuth();
//    
//    // We also have to define certain parameters of the Monte Carlo simulation
//    int num_debris_runs_per_explosion = 1;
//    
////	double altitude_we_explode_until_km = 70;
////    double altitude_we_explode_until_km = alt_explode_until_in;
//    double debris_delta_t = 2.0;
//    
//    // Current Trajectory is probably already in some other mode, find the appropriate StateTimes and switch over to debris mode
//    // Figure out which statetimes we need
//    vector <double> StateTimes;
//    if (StateDot_Option == FIRST_STAGE){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == WHOLE_ROCKET){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == SUBORBITAL) {
//        StateTimes = SuborbitalStateTimes;      }
//    else {
//        cout << "ERROR (MonteCarloDebris): You passed in an unrecognized StateDot_Option.  Exiting...\n\n";
//        exit(2);
//    }
//    
//	int currTraj = 0;
//	double debris_state[8];	//contains state vector at time of breakup PLUS the UTC time
//	
//    // Determine the cutoff_tstep, i.e. when the nominal rocket ascends (for the first time) above altitude_we_explode_until_km
//    std::vector<std::vector<std::vector<double> > > RocketLatLonAltStorage;
//    RocketLatLonAltStorage = TransformToLatLonAlt(1);
//    int rocket_num_time_steps = (INTxx) RocketLatLonAltStorage[0].size();
//    
//    
//    cout << "BIGASS WARNING!!!  Bypassing alt_explode_until in favor of direction specifying the time" << endl;
////	int cutoff_tstep;
////	for (cutoff_tstep = 0; cutoff_tstep < rocket_num_time_steps; cutoff_tstep++) {
////		if (RocketLatLonAltStorage[currTraj][cutoff_tstep][2] > altitude_we_explode_until_km) {
////			break; } }
////	cutoff_tstep--;		//back it off one
//    
//    int cutoff_tstep;
//    cutoff_tstep = ceil(time_explode_until_in/ delta_t );
//    
//    // Switch over
//    StateDot_Option = DEBRIS;
//	Debris TestCatalog;
//    
//    int backwardsStep = 500;
//	// Go backwards so that making footprints is faster
//	for (int tstep = cutoff_tstep; tstep >= 0; tstep = tstep - backwardsStep) {
//		cout << "current debris timestep = " << tstep << endl;
//        cout << "altitude at this timestep = " << RocketLatLonAltStorage[currTraj][tstep][2] << endl;
//        
//		// Copy the state vector into debris_state
//		memcpy((void *) debris_state, (void *) &State_Vector_Storage_Vec[currTraj][tstep][0], 7*sizeof(double));
//		debris_state[7] = Initial_UTC + StateTimes[tstep]/(24*3600.);	// And append the UTC time
//		
//		// Blow things up (generating a new debris profile every time)
//		for (int ix = 0; ix < num_debris_runs_per_explosion; ix++) {
//			TestCatalog.GenerateRandomPieces();
//            //			TestCatalog.ReduceRandomPieces(5.,0.);	//CHOPPING OUT ALL THE HIGH BALLISTIC COEFFICIENTS
//			
//			// Propagate debris
//			initialize_random_atmosphere();	//This isn't really the right thing to do...should have a consistent wind profile for every individual run
//            // What SHOULD be changing here are the debris properties, not the atmosphere.
//            // Though in some ways this gives us maximum variation which is nice.
//            // function sets variable debris_max_time_steps internally
//            
//			DebrisLatLonAltStorage = Propagate_Debris_From_Catalog(debris_state, TestCatalog, debris_delta_t, (unsigned int) tstep);
//            
////            cout << "USING FRISCO DEBRIS!!!!!  HIGHLY EXPERIMENTAL!!!!!  WATCH OUT AND DONT FORGET ABOUT THIS!!!!!~~~~~~~~~~~~~~~~~" << endl;
////			DebrisLatLonAltStorage = Frisco_Propagate_Debris_From_Catalog(debris_state, TestCatalog, debris_delta_t, (unsigned int) tstep);
//            
//            // Debugging here
//            char TestFile[] = "GeneratedFiles/TestDebrisFall.kml";
//            int printThisMany = 20;
//            write_to_google_earth_native(TestFile, printThisMany);
//			// End Debug
//            
//			// Update debris variables
//			DebrisInitialUTC = debris_state[7];
//            
//            // MOVING THIS INSIDE Propagate_Debris_From_Catalog
//            //			//I know that I'm using a constant time step so I can generate DebrisStateTimes on my own
//            //			if (DebrisStateTimes == 0) {    // Is this the right way to do that?  I doubt it.
//            //				delete DebrisStateTimes; }
//            //
//            //            // DebrisStateTimes gets used within function that assembles the total_points_at vector
//            //			DebrisStateTimes = new double [debris_max_time_steps];
//            //			for (int debrisDT = 0; debrisDT < debris_max_time_steps; debrisDT++) {
//            //				DebrisStateTimes[debrisDT] = debrisDT*debris_delta_t; }
//            ////			DebrisFinalUTC = DebrisInitialUTC + DebrisStateTimes[debris_max_time_steps-1]/(3600.*24);
//            
//			vector<vector<Point> > total_points_at = assemble_all_points_debris(tstepMinutes);
//            
//            //            my_footprint.GridTheSky(total_points_at);
//            
//            
//            cout << "COMMENTED OUT THE PART THAT APPENDS DEBRIS TO EXISTING FOOTPRINTS!!!!!  WATCH OUT AND DONT FORGET ABOUT THIS!!!!!~~~~~~~~~~~~~~~~~" << endl;
////            my_footprint.append_to_existing_footprint(total_points_at, DebrisInitialUTC, tstepMinutes, launchLat, launchLon, launchAzimuth);
//            
//		} }
//    
//    return;
//    
//}
//
//
//
//
//
//
//
//
//
//
//
//
///*! MonteCarloDebris: Currently, this function
// *  * Inputs a pre-existing footprint object and updates it
// *  *  * Might be a nice idea to have option allowing the footprint object to be updated
// *  *  * Would also be nice to pass in the state-vectors you want to integrate, rather than figuring them out in this function
// *  * propagates debris to the ground starting from altitude_we_explode_until_km at every tstep = tstep - backwardsStep
// *  * at each considered timestep (of initial explosion) APPENDS the generated points into the existing envelope
// *  * when finished, stores the footprint as a vector in a binary file via store_footprint_as_vector
// */
//void Trajectory::MonteCarloDebris(SkyGrid &my_SkyGrid, double time_explode_until_in){
//    // By this point, the parameters of the architecture should be already defined.  Needed params are:
////    double binSizeKm = 5.;  //Actually, everything in Trajectory should be agnostic to this one
//    
//    // These will be used in the merge steps, passed along
//    double tstepMinutes = my_SkyGrid.getDeltaT();
//    double launchLat = my_SkyGrid.getLaunchLat();
//    double launchLon = my_SkyGrid.getLaunchLon();
//    double launchAzimuth = my_SkyGrid.getLaunchAzimuth();
//    
//    // We also have to define certain parameters of the Monte Carlo simulation
//    int num_debris_runs_per_explosion = 1;
//	double altitude_we_explode_until_km = 100;   //was 50
//    double debris_delta_t = 2.0;
//
//    // Current Trajectory is probably already in some other mode, find the appropriate StateTimes and switch over to debris mode
//    // Figure out which statetimes we need
//    vector <double> StateTimes;
//    if (StateDot_Option == FIRST_STAGE){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == WHOLE_ROCKET){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == SUBORBITAL) {
//        StateTimes = SuborbitalStateTimes;      }
//    else {
//        cout << "ERROR (MonteCarloDebris): You passed in an unrecognized StateDot_Option.  Exiting...\n\n";
//        exit(2);
//    }
//    
//	int currTraj = 0;
//	double debris_state[8];	//contains state vector at time of breakup PLUS the UTC time
//	
//    // Determine the cutoff_tstep, i.e. when the nominal rocket ascends (for the first time) above altitude_we_explode_until_km
//    std::vector<std::vector<std::vector<double> > > RocketLatLonAltStorage;
//    RocketLatLonAltStorage = TransformToLatLonAlt(1);
//    int rocket_num_time_steps = (INTxx) RocketLatLonAltStorage[0].size();
//    
////	int cutoff_tstep;
////	for (cutoff_tstep = 0; cutoff_tstep < rocket_num_time_steps; cutoff_tstep++) {
////		if (RocketLatLonAltStorage[currTraj][cutoff_tstep][2] > altitude_we_explode_until_km) {
////			break; } }
////	cutoff_tstep--;		//back it off one
//    
//    int cutoff_tstep;
//    cutoff_tstep = ceil(time_explode_until_in/ delta_t );
//
//    
////	cout << "DEBUGGING CUTOFF_TSTEP WITH HARDCODE!!!!!!\n";
////	cutoff_tstep = 70;
//    
//    // Switch over
//    StateDot_Option = DEBRIS;
//	Debris TestCatalog;
//    
//    int backwardsStep = 1000;
//	// Go backwards so that making footprints is faster
//	for (int tstep = cutoff_tstep; tstep >= 0; tstep = tstep - backwardsStep) {
//		cout << "current debris timestep = " << tstep << endl;
//        
//		// Copy the state vector into debris_state
//		memcpy((void *) debris_state, (void *) &State_Vector_Storage_Vec[currTraj][tstep][0], 7*sizeof(double));
//		debris_state[7] = Initial_UTC + StateTimes[tstep]/(24*3600.);	// And append the UTC time
//		
//		// Blow things up (generating a new debris profile every time)
//		for (int ix = 0; ix < num_debris_runs_per_explosion; ix++) {
//			TestCatalog.GenerateRandomPieces();
//            //			TestCatalog.ReduceRandomPieces(5.,0.);	//CHOPPING OUT ALL THE HIGH BALLISTIC COEFFICIENTS
//			
//			// Propagate debris
//			initialize_random_atmosphere();	//This isn't really the right thing to do...should have a consistent wind profile for every individual run
//            // What SHOULD be changing here are the debris properties, not the atmosphere.
//            // Though in some ways this gives us maximum variation which is nice.
//            // function sets variable debris_max_time_steps internally
//            
//			DebrisLatLonAltStorage = Propagate_Debris_From_Catalog(debris_state, TestCatalog, debris_delta_t, (unsigned int) tstep);
//            
//			// Update debris variables
//			DebrisInitialUTC = debris_state[7];
//            
//            // MOVING THIS INSIDE Propagate_Debris_From_Catalog
////			//I know that I'm using a constant time step so I can generate DebrisStateTimes on my own
////			if (DebrisStateTimes == 0) {    // Is this the right way to do that?  I doubt it.
////				delete DebrisStateTimes; }
////            
////            // DebrisStateTimes gets used within function that assembles the total_points_at vector
////			DebrisStateTimes = new double [debris_max_time_steps];
////			for (int debrisDT = 0; debrisDT < debris_max_time_steps; debrisDT++) {
////				DebrisStateTimes[debrisDT] = debrisDT*debris_delta_t; }
//////			DebrisFinalUTC = DebrisInitialUTC + DebrisStateTimes[debris_max_time_steps-1]/(3600.*24);
//            
//			vector<vector<Point> > total_points_at = assemble_all_points_debris(tstepMinutes);
//            
////            my_footprint.GridTheSky(total_points_at);
//
////            my_SkyGrid.appendToExistingSkyGrid(total_points_at, DebrisInitialUTC);
////            my_SkyGrid.incorporateNewPoints(total_points_at, DebrisInitialUTC, tstepMinutes, launchLat, launchLon, launchAzimuth);
//            
//            cout << "COMMENTED OUT THE PART THAT APPENDS DEBRIS TO EXISTING FOOTPRINTS!!!!!  WATCH OUT AND DONT FORGET ABOUT THIS!!!!!~~~~~~~~~~~~~~~~~" << endl;
////			my_footprint.append_to_existing_footprint(total_points_at, DebrisInitialUTC, tstepMinutes, launchLat, launchLon, launchAzimuth);
//            
//		} }
//
//    return;
//
//}
//
//
//
//
//
//
//
//
//
//vector<vector<Point> > Trajectory::assemble_all_points_debris(double tstepMinutes) {
//	// Writes the points assuming instantaneous.
//	// In the future, could have a lookahead time and group debris locations with that (note: that's not binning).
//	
//	// Dump those points into the all_points vector (initialized to accomodate all (desired) time steps
//	//Collect all the points
//	//Allocate the vector of points that we're going to write to file
//    int debris_max_time_steps = (INTxx) DebrisStateTimes.size();
//	double timeInFlight = (DebrisStateTimes[debris_max_time_steps-1] - DebrisStateTimes[0])/60.;	//minutes
//	int time_steps_out = (INTxx) ceil(timeInFlight/tstepMinutes);
//	
//	std::vector<std::vector<Point> > total_points_at;
//	total_points_at.assign(time_steps_out,std::vector<Point>());
//	
//	
//    //	// Open a file to write a ground signature for debugging / validation purposes
//    //	ofstream gsig;
//    //	gsig.open("/Users/marian/Documents/MATLAB/tjc_ground_debug2.m",ios::out);
//    //	ofstream gxyz;
//    //	gxyz.open("/Users/marian/Documents/MATLAB/tjc_xyz_debug.m",ios::out);
//    //	cout << "HEY!!!  You're secretly writing a ground signature to a file.  Delete this when done debugging~~~~~~~~~~~~~~" << endl;
//    //	gsig << "tjc_ground = [";
//    //	gxyz << "tjc_xyz = [";
//	
//	
//	int numPieces = (INTxx) DebrisLatLonAltStorage.size();
//	for (int pc = 0; pc < numPieces; pc++) {
//		int numCurrTimeSteps = (INTxx) DebrisLatLonAltStorage[pc].size();
//		for (int tstep = 0; tstep < numCurrTimeSteps; tstep++) {
//			// Figure out the timestep index we're going to bin this into
//			int binTimeStep = (INTxx) floor(DebrisStateTimes[tstep]/(60.*tstepMinutes));
//			
//			// I recall the lon == 0 condition fixed some problem when not all points were initialized.
//			//		I don't know if this is still a problem or not.  Keep it for now.
//			if (!(DebrisLatLonAltStorage[pc][tstep][1] == 0)
//				&& (binTimeStep < time_steps_out)) {
//                double gdlat = DebrisLatLonAltStorage[pc][tstep][0];
//                double lon = DebrisLatLonAltStorage[pc][tstep][1];
//                double curAlt = DebrisLatLonAltStorage[pc][tstep][2];
//                double temp_local_R = get_local_earth_radius(gdlat);
//                Point temp_pt(gdlat, lon, curAlt, temp_local_R);    //TJCHERE
//				
////				Point temp_pt;
////				double temp_local_R = get_local_earth_radius(DebrisLatLonAltStorage[pc][tstep][0]);
////				
////				temp_pt.set_x(DebrisLatLonAltStorage[pc][tstep][1]*R_equator);
////				temp_pt.set_y(DebrisLatLonAltStorage[pc][tstep][0]*temp_local_R);
////				temp_pt.set_z(DebrisLatLonAltStorage[pc][tstep][2]);
////				temp_pt.set_R_local(temp_local_R);
//				
//				total_points_at[binTimeStep].push_back(temp_pt); }
//            
//            //			if (tstep == numCurrTimeSteps-1) {
//            //				gsig << DebrisLatLonAltStorage[pc][tstep][0]*180/PI << "  "
//            //					<< DebrisLatLonAltStorage[pc][tstep][1]*180/PI << "  "
//            //					<< DebrisLatLonAltStorage[pc][tstep][2]*1e3  << endl;
//            //
//            //				gxyz << total_points_at[binTimeStep].back().get_x() << "  "
//            //				<< total_points_at[binTimeStep].back().get_y() << "  "
//            //				<< total_points_at[binTimeStep].back().get_z() << "  "
//            //					<< endl;
//            //			}
//            
//		} }
//    //	gsig << "];\n\n\n";
//    //	gsig.close();
//    //
//    //	gxyz << "];\n\n\n";
//    //	gxyz.close();
//	
//    
////	if (outFileName != "none") {
////		// Now write the stuff to a file
////		ofstream outfile;
////		outfile.open(outFileName.c_str(),ios::out | ios::binary);
////		
////		// Write timing info
////		outfile.write((char *) &time_steps_out, sizeof(time_steps_out));
////		outfile.write((char *) &DebrisInitialUTC, sizeof(DebrisInitialUTC));
////		outfile.write((char *) &tstepMinutes, sizeof(tstepMinutes));
////		
////		// Write vector structure info
////		int num_points_here;
////		for (int t = 0; t < time_steps_out; t++) {
////			num_points_here = (INTxx) total_points_at[t].size();
////			outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
////		
////		// Write the points
////		for (int t = 0; t < time_steps_out; t++) {
////			outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
////		
////		outfile.close();
////	}
//	
//	return total_points_at;
//}
//














// NOTE: Possibly useful but not using as of June 2013 -- LACKS LAUNCH SITE DATA!!!!!
//void Trajectory::write_trajectory_to_file(char *outputFileName, int num_to_write){
//	//Need to be careful when writing to a binary.  the elements of a std::vector are allocated such that they
//	//  are memory contiguous, but the contiguity of a vector of a vector is NOT guaranteed.
//	
////	int num_to_write = 1;	//would be num_per_batch if we were writing a bunch
//
////    std::vector<std::vector<std::vector<double> > > LatLonAltStorage = TransformToLatLonAlt(num_to_write);
//	
//	// Now write the latlonalt coords to a binary file
//	ofstream outfile;
////	outfile.open("GeneratedFiles/trajectory_points.dat",ios::out | ios::binary);
//	outfile.open(outputFileName,ios::out | ios::binary);
//	
//	// Write the number of distinct trajectories to expect
//	outfile.write((char *) &num_to_write, sizeof(num_to_write));
//	
//    // For each trajectory, write the number of time steps to expect
//	int num_points_here;
//	for (int ix = 0; ix < num_to_write; ix++) {
//		num_points_here = (INTxx) State_Vector_Storage_Vec[ix].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here));  }
//    
//    
////	// For each trajectory, write the number of time steps to expect
////	int num_points_here;
////	for (int ix = 0; ix < num_to_write; ix++) {
////		num_points_here = (INTxx) LatLonAltStorage[ix].size();
////		outfile.write((char *) &num_points_here, sizeof(num_points_here));
////	}
//	
////	// Write the LatLonAlt for each trajectory at every available timestep
////	for (int ix = 0; ix < num_to_write; ix++) {
////		num_points_here = (INTxx) LatLonAltStorage[ix].size();
////		for (int jx = 0; jx < num_points_here; jx++) {
////			for (int kx = 0; kx < 3; kx++) {
////				outfile.write((char *) &LatLonAltStorage[ix][jx][kx], sizeof(double)); } } } 
//	
//	// Now write the basic timing information
//	outfile.write((char *) &Initial_UTC, sizeof(Initial_UTC));
//    
//    int StateTimesLength = (INTxx) RocketStateTimes.size();
//	outfile.write((char *) &StateTimesLength, sizeof(StateTimesLength));
//	
//	// Write the StateTimes array (simple array of doubles, i expect this to be memory contiguous)
//	outfile.write((char *) &RocketStateTimes[0], StateTimesLength*sizeof(double));
//	
//	// Write the size of the state vector
//	int temp_size = state_vec_size;		//Compiler complains about passing reference of static const
//	outfile.write((char *) &temp_size, sizeof(temp_size));
//
//	// Now write the state vectors in ECI at every timestep
//	for (int ix = 0; ix < num_to_write; ix++) {
//		num_points_here = (INTxx) State_Vector_Storage_Vec[ix].size();
//		for (int jx = 0; jx < num_points_here; jx++) {
//			for (int kx = 0; kx < state_vec_size; kx++) {
//				outfile.write((char *) &State_Vector_Storage_Vec[ix][jx][kx], sizeof(double)); } } }
//    
//    // Write the falling stage state times
//    StateTimesLength = (INTxx) StageDownStateTimes.size();
//	outfile.write((char *) &StateTimesLength, sizeof(StateTimesLength));
//    
//    // For each trajectory, write the number of time steps to expect
//	for (int ix = 0; ix < num_to_write; ix++) {
//		num_points_here = (INTxx) FallingStageStateVectorStorage[ix].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here));  }
//
//    // Write the falling stage state vectors
//	for (int ix = 0; ix < num_to_write; ix++) {
//		num_points_here = (INTxx) FallingStageStateVectorStorage[ix].size();
//		for (int jx = 0; jx < num_points_here; jx++) {
//			for (int kx = 0; kx < state_vec_size; kx++) {
//				outfile.write((char *) &FallingStageStateVectorStorage[ix][jx][kx], sizeof(double)); } } }
//    
//    // Write the structural mass jetissoned as the first stage burns out
//	outfile.write((char *) &StageStructMass[0], sizeof(double));
//    
//	outfile.close();
//    
//    //TODO: Generalize the stage mass stuff to accomodate multiple stages
//    //TODO: Write the individual stage time lengths
//	
//	return;
//}

// NOTE: Possibly useful but not using as of June 2013 -- LACKS LAUNCH SITE DATA!!!!!
//void Trajectory::read_trajectory_from_file(char *outputFileName){
//	//Need to be careful when writing to a binary.  the elements of a std::vector are allocated such that they
//	//  are memory contiguous, but the contiguity of a vector of a vector is NOT guaranteed.
//	
//	//	int num_to_write = 1;	//would be num_per_batch if we were writing a bunch
//	
//	ifstream infile;
//	infile.open(outputFileName,ios::in | ios::binary);
//	
//	std::vector<double> ZeroVec3Size(3,0.);
//	
//	int num_trajectories;
//	if (infile.good()) {
//		
//		infile.read((char *) &num_trajectories, sizeof(num_trajectories));	// read in how many time steps there are
//        
////		RocketLatLonAltStorage.assign(num_trajectories,std::vector<std::vector<double> >() );
//        
//		cout << "num_trajectories = " << num_trajectories << endl;
//		State_Vector_Storage_Vec.assign(num_trajectories,std::vector<std::vector<double> >() );
//        FallingStageStateVectorStorage.assign(num_trajectories,std::vector<std::vector<double> >() );
//        
//		int num_points_here;
//        double pointsHereVec [num_trajectories];
//		for (int ix = 0; ix < num_trajectories; ix++) {
//			infile.read((char *) &num_points_here, sizeof(num_points_here));
//			cout << "num_point_here = " << num_points_here << endl;
//            pointsHereVec[ix] = num_points_here; }
////			State_Vector_Storage_Vec[ix].assign(num_points_here,ZeroVec3Size); }
//        
//        
//        
//        
////		for (int ix = 0; ix < num_trajectories; ix++) {
////			infile.read((char *) &num_points_here, sizeof(num_points_here));
////			cout << "num_point_here = " << num_points_here << endl;
////			RocketLatLonAltStorage[ix].assign(num_points_here,ZeroVec3Size);}
//        
////		//Read in each value and assign it
////		for (int ix = 0; ix < num_trajectories; ix++) {
////			num_points_here = RocketLatLonAltStorage[ix].size();
////			cout << "num_point_here = " << num_points_here << endl;
////			for (int jx = 0; jx < num_points_here; jx++) {
////				for (int kx = 0; kx < 3; kx++) {
////					infile.read((char *) &RocketLatLonAltStorage[ix][jx][kx], sizeof(double)); } } }
//		
//		// Now read the basic timing information
//		infile.read((char *) &Initial_UTC, sizeof(Initial_UTC));
//        
//        int num_time_steps;
//		infile.read((char *) &num_time_steps, sizeof(num_time_steps));
//		
//		// Allocate the StateTimes array (would rather have this be a std::vector...try that later)
////		RocketStateTimes = new double[num_time_steps];
//        RocketStateTimes.assign(num_time_steps,0.);
//		
//		// Read the StateTimes array (simple array of doubles, i expect this to be memory contiguous)
//		infile.read((char *) &RocketStateTimes[0], num_time_steps*sizeof(double));
//        for (int ix = 0; ix < num_time_steps; ix++){
//            cout << "rocetstatetimes = " << RocketStateTimes[ix] << endl;
//        }
//		
//		// Read the size of the state vector This is actually kinda pointless
//        int temp_vec_size;
//		infile.read((char *) &temp_vec_size, sizeof(temp_vec_size));
//		std::vector<double> ZeroVecStateSize(state_vec_size,0.);
//        
//		// Now read the state vectors in ECI at every timestep
//		State_Vector_Storage_Vec.assign(num_trajectories,std::vector<std::vector<double> >());
//		for (int ix = 0; ix < num_trajectories; ix++) {
////			num_points_here = RocketLatLonAltStorage[ix].size();
//			State_Vector_Storage_Vec[ix].assign(pointsHereVec[ix], ZeroVecStateSize);
//			for (int jx = 0; jx < pointsHereVec[ix]; jx++) {
//				for (int kx = 0; kx < state_vec_size; kx++) {
//					infile.read((char *) &State_Vector_Storage_Vec[ix][jx][kx], sizeof(double)); } } }
//        
//        
////        // Write the falling stage state times
////        StateTimesLength = (INTxx) StageDownStateTimes.size();
////        outfile.write((char *) &StateTimesLength, sizeof(StateTimesLength));
//        infile.read((char *) &num_time_steps, sizeof(num_time_steps));
////        
////        // For each trajectory, write the number of time steps to expect
////        for (int ix = 0; ix < num_to_write; ix++) {
////            num_points_here = (INTxx) FallingStageStateVectorStorage[ix].size();
////            outfile.write((char *) &num_points_here, sizeof(num_points_here));  }
//        for (int ix = 0; ix < num_trajectories; ix++) {
//            infile.read((char *) &num_points_here, sizeof(num_points_here));
//            FallingStageStateVectorStorage[ix].assign(num_points_here, ZeroVecStateSize);  }
////
////        // Write the falling stage state vectors
////        for (int ix = 0; ix < num_to_write; ix++) {
////            num_points_here = (INTxx) FallingStageStateVectorStorage[ix].size();
////            for (int jx = 0; jx < num_points_here; jx++) {
////                for (int kx = 0; kx < state_vec_size; kx++) {
////                    outfile.write((char *) &FallingStageStateVectorStorage[ix][jx][kx], sizeof(double)); } } }
//        
//        for (int ix = 0; ix < num_trajectories; ix++) {
//            num_points_here = (INTxx) FallingStageStateVectorStorage[ix].size();
//            for (int jx = 0; jx < num_points_here; jx++) {
//                for (int kx = 0; kx < state_vec_size; kx++) {
//                    infile.read((char *) &FallingStageStateVectorStorage[ix][jx][kx], sizeof(double));  } } }
//
//        StageStructMass = new double [1];
//		infile.read((char *) &StageStructMass[0], sizeof(double));
//		
//        cout << "first_state_struct-mas = " << StageStructMass[0] << endl;
//        
//        //		cout << "reading in Storage[0][0][0] = " << RocketLatLonAltStorage[0][0][0]*180/PI << endl;
//        
//		
//	} else {
//		cout << "Opening trajectory binary file failed" << endl;
//	}
//	
//    //	for (int ix = 0; ix < RocketLatLonAltStorage[0].size(); ix++) {
//    //		cout << "ix = " << ix << "   gdlat = " << RocketLatLonAltStorage[0][ix][0] << endl; }
//	
//	
//	infile.close();
//	
//	
//	
//	
//	return;
//}



//// NOTE: I don't seem to be using this at the moment June 2013 but could be useful -- LACKS LAUNCH INFO
//void Trajectory::write_debris_to_file(unsigned int idNum){
//	//Need to be careful when writing to a binary.  the elements of a std::vector are allocated such that they
//	//  are memory contiguous, but the contiguity of a vector of a vector is NOT guaranteed.
//	
//    int debris_max_time_steps = (INTxx) DebrisStateTimes.size();
//
//	// Convert everything to LatLon
//	// Find the transformation matrices ECEF__C__ECI at every possible timestep	
//	MatrixM Zero3x3 (3,3,0.);
//	std::vector< MatrixM > ECEF__C__ECI(debris_max_time_steps,Zero3x3);
//	for (int i = 0; i < debris_max_time_steps; i++) {
//		ECEF__C__ECI[i] = trans (ECEF_to_ECI_rotation(Initial_UTC + DebrisStateTimes[i]/(24.*3600)) ); }
//	
//	
//	//First generate the vector of vectors as we will want them to be read in
//	ColVec Geodetic(3,0.);
//	ColVec ECI(3,0.);
//	std::vector<double> ZeroVec3Size(3,0.);
//	std::vector<std::vector<std::vector<double> > > LatLonAltStorage(num_debris_pieces,std::vector<std::vector<double> >() );
//	// Load the geodetic coordinates into the overall storage vector
//	for (int j = 0; j < num_debris_pieces; j++) {
//		
//		int i_limit = (INTxx) State_Vector_Storage_Vec[j].size();
//		LatLonAltStorage[j].assign(i_limit,ZeroVec3Size);
//		
//		for (int i = 0; i < i_limit; i++) {
//			ECI(0) = State_Vector_Storage_Vec[j][i][0];
//			ECI(1) = State_Vector_Storage_Vec[j][i][1];
//			ECI(2) = State_Vector_Storage_Vec[j][i][2];
//			
//			Geodetic = ECEF_To_Geodetic( prod(ECEF__C__ECI[i],ECI));
//			
//			LatLonAltStorage[j][i][0] = Geodetic(0);	//radians
//			LatLonAltStorage[j][i][1] = Geodetic(1);	//radians
//			LatLonAltStorage[j][i][2] = Geodetic(2);	//km
//		} }
//	
//	
//	
//	// Now write the latlonalt coords to a binary file
//	string fileString = "GeneratedFiles/debris_points";
//	string extString = ".dat";
//	char buffer[5];
//	sprintf(buffer,"%i",idNum);
//	fileString = fileString + buffer + extString;
//	
//	ofstream outfile;
//	outfile.open(fileString.c_str(),ios::out | ios::binary);
//	
//	outfile.write((char *) &num_debris_pieces, sizeof(num_debris_pieces));
//	
//	int num_points_here;
//	for (int ix = 0; ix < num_debris_pieces; ix++) {
//		num_points_here = (INTxx) LatLonAltStorage[ix].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here));
//	}
//		
//	//int vector_size = LatLonAltStorage.capacity() + num_debris_pieces*LatLonAltStorage[0].capacity();	
//	//cout << "vector size is probably = " << vector_size << endl;
//
//	for (int ix = 0; ix < num_debris_pieces; ix++) {
//		num_points_here = (INTxx) LatLonAltStorage[ix].size();
//		for (int jx = 0; jx < num_points_here; jx++) {
//			for (int kx = 0; kx < 3; kx++) {
//				outfile.write((char *) &LatLonAltStorage[ix][jx][kx], sizeof(double)); } } } 
//
//	// Now write the basic timing information
//	outfile.write((char *) &Initial_UTC, sizeof(Initial_UTC));
//	outfile.write((char *) &debris_max_time_steps, sizeof(debris_max_time_steps));
//	
//	// Write the StateTimes array (simple array of doubles, i expect this to be memory contiguous)
//	outfile.write((char *) &DebrisStateTimes[0], debris_max_time_steps*sizeof(double));
//	
//	outfile.close();
//		
//	
//	
//	return;
//}



//
//// This function appears to get used to shift a thrust profile from one launch location to another
////  e.g. turning the launch from the cape into a launch from Wallops (MARS)
//void Trajectory::Transform_Nominal(string nominalFile, string outputFile) {
//	Get_Nominal_Trajectory(nominalFile);
//
//	// Coords we want to transform to
//	double MarsLong = -75.4882;
//	double MarsLat = 37.8338;
//    //	double MarsAlt = 0;
//
//	// First transform thrust so that z-axis aligns
//		//Actually, we're in ECI so z's should already be aligned
//
//	// Instead of rotating coord systems, we want to rotate the thrust vector itself, so take transposes
//	// Rotate over in longitude
//	double deltaLong = MarsLong - lon;	//Positive for MARS
//	MatrixM Rot1 = trans(Rotation(3,deltaLong));
//
//	// Rotate up in lat
//	double deltaLat = MarsLat - gdlat;
//	MatrixM Rot2 = trans(Rotation(2,-deltaLat));
//
//	MatrixM FullRotation = prod(Rot2,Rot1);
//
//	ColVec OldThrust(3,0.);
//	ColVec NewThrust(3,0.);
//
//	// Set the precision and make sure output is set to scientific
//	cout.precision(9);
//	cout << std::scientific;
//
//	for (int ix = 0; ix < NumSteps; ix++){
//		OldThrust(0) = ThrustUx[ix];
//		OldThrust(1) = ThrustUy[ix];
//		OldThrust(2) = ThrustUz[ix];
//
//		// Now can rotate the vector
//			//FullRotation
//		NewThrust = prod(FullRotation,OldThrust);
//
//		cout << std::setw(20) << NewThrust(0) << std::setw(20) << NewThrust(1) << std::setw(20) << NewThrust(2) << std::setw(20) << ThrustMag[ix] << std::setw(20) << ThrustTime[ix] << "\n";
//	}
//
//	exit(8);	//Don't let the program continue of you ever wind up in here.  Kill it to be safe.
//	return;
//}












// //This function is out for now, BUT DO NOT GRAVEYARD.
//void Trajectory::write_debris_to_file(unsigned int idNum, std::vector<std::vector<std::vector<double> > > &LatLonAltStorage){
//	//Need to be careful when writing to a binary.  the elements of a std::vector are allocated such that they
//	//  are memory contiguous, but the contiguity of a vector of a vector is NOT guaranteed.
//
//    cout << "If you're using this BE CAREFUL!!!!!!!  PLEASE READ THE COMMENTS FOR THIS FUNCTION!!!" << end;
//    // I'm commenting this file out because I think it's not useful at the moment, however it will become useful again in the future,
//    //  after some modifications.  Should be passing the state times vector in to be explicit.  In fact, maybe should pass everything in.
//    //  Also obviously want to not hard-code the name and location of the points file.
//
//	// Now write the latlonalt coords to a binary file
//	string fileString = "GeneratedFiles/debris_points";
//	string extString = ".dat";
//	char buffer[5];
//	sprintf(buffer,"%i",idNum);
//	fileString = fileString + buffer + extString;
//
//	ofstream outfile;
//	outfile.open(fileString.c_str(),ios::out | ios::binary);
//
//	outfile.write((char *) &num_debris_pieces, sizeof(num_debris_pieces));
//
//	int num_points_here;
//	for (int ix = 0; ix < num_debris_pieces; ix++) {
//		num_points_here = (INTxx) LatLonAltStorage[ix].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here));
//	}
//
//	//int vector_size = LatLonAltStorage.capacity() + num_debris_pieces*LatLonAltStorage[0].capacity();
//	//cout << "vector size is probably = " << vector_size << endl;
//
//	for (int ix = 0; ix < num_debris_pieces; ix++) {
//		num_points_here = (INTxx) LatLonAltStorage[ix].size();
//		for (int jx = 0; jx < num_points_here; jx++) {
//			for (int kx = 0; kx < 3; kx++) {
//				outfile.write((char *) &LatLonAltStorage[ix][jx][kx], sizeof(double)); } } }
//
//	// Now write the basic timing information
//	outfile.write((char *) &Initial_UTC, sizeof(Initial_UTC));
//	outfile.write((char *) &debris_max_time_steps, sizeof(debris_max_time_steps));
//
//	// Write the StateTimes array (simple array of doubles, i expect this to be memory contiguous)
//	outfile.write((char *) &DebrisStateTimes[0], debris_max_time_steps*sizeof(double));
//
//	outfile.close();
//
//	return;
//}










// ~~~~~~~~~~~~~~~~~~~~~ !!!!GRAVEYARD!!!! ~~~~~~~~~~~~~~~~~~~~~~~~









//int Trajectory::write_batch_to_google_earth(string basename, int batch_num) {
//	
//	// Find the transformation matrices ECEF__C__ECI at every possible timestep
//	MatrixM Zero3x3 (3,3,0.);
//	std::vector< MatrixM > ECEF__C__ECI(num_time_steps,Zero3x3);
//	for (int i = 0; i < num_time_steps; i++) {
//		ECEF__C__ECI[i] = trans (ECEF_to_ECI_rotation(Initial_UTC + StateTimes[i]/(24.*3600)) ); }
//	
//	ofstream outfile;
//	
//	
//	
//    //	opt_traj.open(nominal_file.c_str(), ifstream::in);
//    
//    //	basename += "_" + num2str(batch_num);
//	
//	string testbase("/Users/marian/Documents/Research/GoogleEarth/Trajectory");
//    
//    //	string testbase("/Users/marian/Documents/Research/GoogleEarth/Trajectory");
//	char buffer [5];
//	string testfile;
//	
//	ColVec Geodetic(3,0.);
//	ColVec ECI(3,0.);
//	
//	int j_limit = -1;
//	if (StateDot_Option == FIRST_STAGE) {
//		j_limit = num_per_batch; }
//	else if (StateDot_Option == DEBRIS) {
//		j_limit = num_debris_pieces; }
//	else if (StateDot_Option == SUBORBITAL) {
//		j_limit = num_per_batch; }
//    
//	// Write the geodetic coordinates of the rocket trajectories
//	for (int j = 0; j < j_limit; j++) {
//        
//		sprintf(buffer, "%i", j+1);
//        //		testfile = testbase + buffer;
//		testfile = testbase + buffer;
//		cout << "landed here" << endl;
//		
//		outfile.open(testfile.c_str());
//        
//		int i_limit = State_Vector_Storage_Vec[j].size();
//		for (int i = 0; i < i_limit; i++) {
//			
//			ECI(0) = State_Vector_Storage_Vec[j][i][0];
//			ECI(1) = State_Vector_Storage_Vec[j][i][1];
//			ECI(2) = State_Vector_Storage_Vec[j][i][2];
//			
//			Geodetic = ECEF_To_Geodetic( prod(ECEF__C__ECI[i],ECI));
//			
//			outfile << Geodetic(1)*180./PI << "  " << Geodetic(0)*180./PI << "  " << (Geodetic(2)-local_elevation)*1e3 << endl;
//		}
//		
//		outfile.close();
//	}
//	
//	//	outfile << "];\n\n" << "plot3(y(:,1),y(:,2),y(:,3))\n" << "axis equal" << endl;
//	
//	return 1;
//}




//void Trajectory::write_LRHC_points_BlockOn(char *outFileName, double tstepMin, double startUTC, double minutesOn) {
//    // This function writes the LRHC ellipse around the nominal trajectory to a points file so that it may be wrapped up into a footprint
//    // Here, we bin the points according to time
//    //  ASSUMES: Traditional rocket that goes up and never re-enters NAS; will break if you use this for something that re-enters
//    //  ASSUMES: Your first timestep has the rocket within the NAS
//
//	// Do some error checking
//	ofstream outfile;
//	outfile.open(outFileName,ios::out | ios::binary);
//
//	if (outfile.bad()) {
//		cerr << "LRHC FILE FAILED TO OPEN FOR OUTPUT!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//	// Things we'll need
//	Point temp_pt;
//	std::vector<std::vector<Point> > total_points_at;	//[timestep][points]
//
//	int num_range;
//
//    //Allocate the vector of points that we're going to write to file
//    //		double timeInFlight = (RocketStateTimes[rocket_num_time_steps-1] - RocketStateTimes[0])/60.;	//minutes
//    int time_steps_out = (INTxx) ceil(minutesOn/tstepMin);
//    total_points_at.assign(time_steps_out,std::vector<Point>());
//    num_range = time_steps_out;
//
//    // Collecting the points from a trajectory ellipse / tube
//    int num_points_in_trajectory = (INTxx) Single_Trajectory_Tube.size();
//    if (num_points_in_trajectory > 0) {	 //checking to make sure ellipse has been allocated
//        for (int curTimeIX = 0; curTimeIX < time_steps_out; curTimeIX++) {
//            for (int t = 0; t < num_points_in_trajectory; t++) {
//                for (int ix = 0; ix < Single_Trajectory_Tube[t].size(); ix++) {
//                    total_points_at[curTimeIX].push_back(Single_Trajectory_Tube[t][ix]); } } } }
//    else {
//        cout << "You do not appear to have generated the LRHC boundary~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//
//	// Begin writing to file
//	outfile.write((char *) &num_range, sizeof(num_range));
//	outfile.write((char *) &startUTC, sizeof(startUTC));
//	outfile.write((char *) &tstepMin, sizeof(tstepMin));
//
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
//
//	outfile.close();
//
//	// Used for debugging purposes only
//	bool debug = false;
//	if (debug) {
//		string debugFileName("GeneratedFiles/LRHC_PLAINTEXT_DEBUG.txt");
//		ofstream debugfile;
//		debugfile.open(debugFileName.c_str(), ios::out);
//
//		debugfile << num_range << endl << startUTC << endl << tstepMin << endl;
//
//		int num_points_here;
//		for (int t = 0; t < num_range; t++) {
//			num_points_here = (INTxx) total_points_at[t].size();
//			for (int ix = 0; ix < num_points_here; ix++) {
//				debugfile << total_points_at[t][ix] << endl; }
//			debugfile << "\n\n\n\n\n\n\n\n\n\n\n\n\n"; }
//
//		debugfile.close();
//	}
//	// End of debugging code
//
//
//	return;
//}












// ---- ALL Points files should have this format -----
// int				number_of_time_steps
// double			UTC_at_start
// double			delta_t_in_minutes
// ints				for (number_of_time_steps) { number_of_points_at_this_time_step }
// vector<Point>	for (number_of_time_steps) { vector_of_points_at_this_time_step }
//void Trajectory::write_LRHC_points(char *outFileName, bool isInstantaneous, double tstepMin, double startUTC = -1, double minutesOn = -1) {
//    // This function writes the LRHC ellipse around the nominal trajectory to a points file so that it may be wrapped up into a footprint
//    // Here, we bin the points according to time
//    //  ASSUMES: Traditional rocket that goes up and never re-enters NAS; will break if you use this for something that re-enters
//    //  ASSUMES: Your first timestep has the rocket within the NAS
//
//	// Do some error checking
//	ofstream outfile;
//	outfile.open(outFileName,ios::out | ios::binary);
//
//	if (outfile.bad()) {
//		cerr << "LRHC FILE FAILED TO OPEN FOR OUTPUT!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//	// Things we'll need
//	Point temp_pt;
//	std::vector<std::vector<Point> > total_points_at;	//[timestep][points]
//
//	int num_range;
//	if (isInstantaneous){
//		// Instead of a solid shape that turns on for a period of time, this will be an ellipse (roughly) propagating upward like a smoke ring
//		// Time binning into minute intervals...may actually make this specific option useless :(
//
//		// We'll want to output this, so save it
//		startUTC = Initial_UTC;
//
//		//Allocate the vector of points that we're going to write to file
//		double timeInFlight = (StateTimes[num_time_steps-1] - StateTimes[0])/60.;	//minutes
//		int time_steps_out = (INTxx) ceil(timeInFlight/tstepMin);
//		total_points_at.assign(time_steps_out,std::vector<Point>());
//		num_range = time_steps_out;
//
//		// Collecting the points from a trajectory ellipse / tube
//		int num_points_in_trajectory = (INTxx) Single_Trajectory_Tube.size();
//		if (num_points_in_trajectory > 0) {										//checking to make sure ellipse has been allocated
//			for (int t = 0; t < num_points_in_trajectory; t++) {
//				// NEED SOME LOGIC HERE TO DETERMINE WHAT TIME BIN THE CURRENT RING OF POINTS WILL GO IN
//				int curTimeIX = (INTxx) floor(StateTimes[t]/(60*tstepMin));
//				if (curTimeIX < time_steps_out) {
//					for (int ix = 0; ix < Single_Trajectory_Tube[t].size(); ix++) {
//						total_points_at[curTimeIX].push_back(Single_Trajectory_Tube[t][ix]); } }
//				else {
//					cerr << "SOMETHING IS WRONG IN TIME-BINNING THE LRHC ELLIPSE!!!~~~~~~~~~~~~~~~~~" << endl;} } }
//		else {
//			cout << "You do not appear to have generated the LRHC boundary~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//	} //ends if(isInstantaneous)
//
//	else {
//		// Going to turn on the full ellipse tube for a specified amount of time, starting at the specified time
//		if ((startUTC == -1) || (minutesOn == -1)) {
//			cerr << "YOU HAVE TO SUPPLY startUTC and minutesOn!!!!!!!!!!~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//			exit(8943234); }
//
//		//Allocate the vector of points that we're going to write to file
//        //		double timeInFlight = (RocketStateTimes[num_time_steps-1] - RocketStateTimes[0])/60.;	//minutes
//		int time_steps_out = (INTxx) ceil(minutesOn/tstepMin);
//		total_points_at.assign(time_steps_out,std::vector<Point>());
//		num_range = time_steps_out;
//
//		// Collecting the points from a trajectory ellipse / tube
//		int num_points_in_trajectory = (INTxx) Single_Trajectory_Tube.size();
//		if (num_points_in_trajectory > 0) {	 //checking to make sure ellipse has been allocated
//			for (int curTimeIX = 0; curTimeIX < time_steps_out; curTimeIX++) {
//				for (int t = 0; t < num_points_in_trajectory; t++) {
//                    for (int ix = 0; ix < Single_Trajectory_Tube[t].size(); ix++) {
//                        total_points_at[curTimeIX].push_back(Single_Trajectory_Tube[t][ix]); } } } }
//		else {
//			cout << "You do not appear to have generated the LRHC boundary~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//	} //ends else isInstantaneous
//
//
//	// Begin writing to file
//	outfile.write((char *) &num_range, sizeof(num_range));
//	outfile.write((char *) &startUTC, sizeof(startUTC));
//	outfile.write((char *) &tstepMin, sizeof(tstepMin));
//
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
//
//	outfile.close();
//
//	// Used for debugging purposes only
//	bool debug = false;
//	if (debug) {
//		string debugFileName("GeneratedFiles/LRHC_PLAINTEXT_DEBUG.txt");
//		ofstream debugfile;
//		debugfile.open(debugFileName.c_str(), ios::out);
//
//		debugfile << num_range << endl << startUTC << endl << tstepMin << endl;
//
//		int num_points_here;
//		for (int t = 0; t < num_range; t++) {
//			num_points_here = (INTxx) total_points_at[t].size();
//			for (int ix = 0; ix < num_points_here; ix++) {
//				debugfile << total_points_at[t][ix] << endl; }
//			debugfile << "\n\n\n\n\n\n\n\n\n\n\n\n\n"; }
//
//		debugfile.close();
//	}
//	// End of debugging code
//
//
//	return;
//}

/*!
 * This function calls MonteCarloDebris in such a way that the original footprint object is left untouched and the resulting
 *  footprint is instead dumped to a file.
 *
 */
//void Trajectory::MonteCarloDebris(Footprint3D &my_footprint, string footprintVectorFile){
//    Footprint3D copyFootprint = my_footprint;
//
//    MonteCarloDebris(copyFootprint);
//    copyFootprint.store_footprint_as_vector(footprintVectorFile);  //must also save the timing information found in points files
//
//    return;
//}






//~~~~~~~~~~~~~~~~~~~~~~ Functions that write data to files ~~~~~~~~~~~~~~~~~~~~~~~~~~~

//int Trajectory::write_batch_to_file() {
//
//    double *StateTimes;
//    if (StateDot_Option == FIRST_STAGE){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == SUBORBITAL) {
//        StateTimes = SuborbitalStateTimes;      }
//
//	ofstream outfile;
//	outfile.open("/Users/marian/Documents/MATLAB/plot_trajectory.m");
//	outfile << "clear all; close all; clc;" << "\n\n";
//	outfile << "DU = " << DU << endl;
//
//	// Write the inertial postions of the rocket trajectories
//	for (int j = 0; j < num_per_batch; j++) {
//		outfile << "O_r_S__eci{" << j+1 << "} = [\n" << endl;
//		for (int i = 0; i < rocket_num_time_steps; i++) {
//			outfile << State_Vector_Storage_Vec[j][i].at(0) << "  " << State_Vector_Storage_Vec[j][i].at(1) << "  " << State_Vector_Storage_Vec[j][i].at(2) << endl;
//		}
//		outfile << "];\n\n";
//	}
//
//	// Find the transformation matrices ECEF__C__ECI at every timestep
//	MatrixM Zero3x3 (3,3,0.);
//	std::vector< MatrixM > ECEF__C__ECI(rocket_num_time_steps,Zero3x3);
//	for (int i = 0; i < rocket_num_time_steps; i++) {
//		ECEF__C__ECI[i] = trans( (ECEF_to_ECI_rotation(Initial_UTC + StateTimes[i]/(24.*3600))) ); }
//
//
//	ColVec Geodetic(3,0.);
//	ColVec ECI(3,0.);
//	// Write the geodetic coordinates of the rocket trajectories
//	for (int j = 0; j < num_per_batch; j++) {
//		outfile << "O_r_S__geodetic{" << j+1 << "} = [\n" << endl;
//		for (int i = 0; i < rocket_num_time_steps; i++) {
//
//			ECI(0) = State_Vector_Storage_Vec[j][i][0];
//			ECI(1) = State_Vector_Storage_Vec[j][i][1];
//			ECI(2) = State_Vector_Storage_Vec[j][i][2];
//
////			cout << "i = " << i << "  j = " << j << "   x = " << State_Vector_Storage_Vec[j][i][0] << endl;
//
//			Geodetic = ECEF_To_Geodetic( prod(ECEF__C__ECI[i],ECI) );
////			cout << "ECI = " << ECI << endl;
//
//			outfile << Geodetic(0) << "  " << Geodetic(1) << "  " << Geodetic(2) << endl;
//		}
//		outfile << "];\n\n";
//	}
//
//	// Write the time vector used
//	outfile << "TimeVec = [\n" << endl;
//	for (int i = 0; i < rocket_num_time_steps; i++) {
//		outfile << StateTimes[i] << endl;
//	}
//	outfile << "];\n\n";
//
//
////	outfile << "Thrust_Interp = [\n" << endl;
////	for (int i = 0; i < last_i; i++) {
////		outfile << Thrust_Force_Storage[i].at(0) << "  " << Thrust_Force_Storage[i].at(1) << "  " << Thrust_Force_Storage[i].at(2) << endl;
////	}
////	outfile << "];\n\n";
////
////	outfile << "Drag_Acc_Prop = [\n" << endl;
////	for (int i = 0; i < last_i; i++) {
////		outfile << Drag_Acc_Storage[i].at(0) << "  " << Drag_Acc_Storage[i].at(1) << "  " << Drag_Acc_Storage[i].at(2) << endl;
////	}
////	outfile << "];\n\n";
////
////	outfile << "Mass_Prop = [\n" << endl;
////	for (int i = 0; i < last_i; i++) {
////		outfile << State_Vector_Storage[i].at(6) << endl;
////	}
////	outfile << "];\n\n";
////
////	outfile << "Density_Profile_Prop = [\n" << endl;
////	for (int cur_alt = 0; cur_alt < 300; cur_alt++) {
////		outfile << cur_alt << "   " << DensityCurve(cur_alt) << endl;
////	}
////	outfile << "];\n\n";
//
//	outfile << "\n\n[xeplot, yeplot, zeplot] = ellipsoid(0.0, 0.0, 0.0, DU, DU, DU, 30);";
//
//	outfile << "\n% Setup earth plotting data";
//	outfile << "\nfigure(1)";
//	outfile << "\nsurface(xeplot, yeplot, zeplot, 'FaceColor', 'blue', 'EdgeColor', 'black');";
//
//	outfile << "\nfor i=1:" << num_per_batch;
//	outfile << "\n\thold on;";
//	outfile << "\n\tplot3(O_r_S__eci{i}(:,1), ...";
//	outfile << "\n\t	  O_r_S__eci{i}(:,2), ...";
//	outfile << "\n\t	  O_r_S__eci{i}(:,3), 'red');";
//	outfile << "\nend";
//
//	outfile << "\nview(3)";
//
//	outfile << "\naxis equal;";
//	outfile << "\ntitle('Orbital Decay Due To Atmospheric Drag (perigee_altitude = 150km)')";
//	outfile << "\nxlabel('x [km]')";
//	outfile << "\nylabel('y [km]')";
//	outfile << "\nzlabel('z [km]')";
//
//
//	//	outfile << "];\n\n" << "plot3(y(:,1),y(:,2),y(:,3))\n" << "axis equal" << endl;
//	outfile.close();
//
//	return 1;
//}






//void Trajectory::write_LRHC_points_file(char *outFileName, double tstepMin) {
//    // This function writes the LRHC ellipse around the nominal trajectory to a points file so that it may be wrapped up into a footprint
//    // Here, we bin the points according to time
//    //  ASSUMES: Traditional rocket that goes up and never re-enters NAS; will break if you use this for something that re-enters
//    //  ASSUMES: Your first timestep has the rocket within the NAS
//
//	// Do some error checking
//	ofstream outfile;
//	outfile.open(outFileName,ios::out | ios::binary);
//
//	if (outfile.bad()) {
//		cerr << "LRHC FILE FAILED TO OPEN FOR OUTPUT!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//	// Things we'll need
//	Point temp_pt;
//	std::vector<std::vector<Point> > total_points_at;	//[timestep][points]
//
//    // Figure out which statetimes we need
//    vector <double> StateTimes;
//    if (StateDot_Option == FIRST_STAGE){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == WHOLE_ROCKET){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == SUBORBITAL) {
//        StateTimes = SuborbitalStateTimes;      }
//
//	int num_range;
//    // Instead of a solid shape that turns on for a period of time, this will be an ellipse (roughly) propagating upward like a smoke ring
//    // Time binning into minute intervals...may actually make this specific option useless :(
//
//    // We'll want to output this, so save it
//    double startUTC = Initial_UTC;
//
//    //Allocate the vector of points that we're going to write to file
////    int num_points_in_trajectory = (INTxx) Single_Trajectory_Tube.size();   //poorly named, time points, time steps
//
//    // This function is going to only consider times up to the first staging event
//    double timeInFlight = (StageTimes[1] - StateTimes[0])/60.;	//minutes
//    int time_steps_out = (INTxx) ceil(timeInFlight/tstepMin);
//    total_points_at.assign(time_steps_out,std::vector<Point>());
//    num_range = time_steps_out;
//
//    int numRuns = num_per_batch;
//    std::vector<std::vector<std::vector<double> > > LatLonAltStorage;
//    LatLonAltStorage = TransformToLatLonAlt(numRuns);
//
//    Point tempPoint;
//
//    double altitudeLimit = 18.288;      // [km] This is the top of the NAS
//
//    for (int runNumber = 0; runNumber < numRuns; runNumber++){
//        int numPtsHere = (INTxx) LatLonAltStorage[runNumber].size();
//        for (int t = 0; t < numPtsHere; t++){
//            double curAlt = LatLonAltStorage[runNumber][t][2];
//            if (curAlt <= altitudeLimit) {
//                int curTimeIX = (INTxx) floor(StateTimes[t]/(60*tstepMin));
//
//                double gdlat = LatLonAltStorage[runNumber][t][0];
//                double rLocal = get_local_earth_radius(gdlat);
//
//                // Anticipate some problems here
//                if (curAlt == 0) {
//                    curAlt = 0.001; //Set at 1m above ground to avoid division by zero
//                }
//                else if (curAlt < 0) {
//                    cout << "ERROR: should not have al0 negative in here, come fix this!!!!~~~~~~~~~" << endl;
//                }
//
//                // Location of nominal trajectory on the ground
//                double curX = R_equator * LatLonAltStorage[runNumber][t][1];  // R_equator * longitude_angle
//                double curY = rLocal * gdlat;                // Local Radius of Earth * gdlat
//
//                tempPoint.set_x(curX);
//                tempPoint.set_y(curY);
//                tempPoint.set_z(curAlt);
//                tempPoint.set_R_local(rLocal);
//
//                total_points_at[curTimeIX].push_back(tempPoint);
//            }
//        } }
//
//
//	// Begin writing to file
//	outfile.write((char *) &num_range, sizeof(num_range));
//	outfile.write((char *) &startUTC, sizeof(startUTC));
//	outfile.write((char *) &tstepMin, sizeof(tstepMin));
//
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
//
//	outfile.close();
//
//	// Used for debugging purposes only
//	bool debug = false;
//	if (debug) {
//		string debugFileName("GeneratedFiles/LRHC_PLAINTEXT_DEBUG.txt");
//		ofstream debugfile;
//		debugfile.open(debugFileName.c_str(), ios::out);
//
//		debugfile << num_range << endl << startUTC << endl << tstepMin << endl;
//
//		int num_points_here;
//		for (int t = 0; t < num_range; t++) {
//			num_points_here = (INTxx) total_points_at[t].size();
//			for (int ix = 0; ix < num_points_here; ix++) {
//				debugfile << total_points_at[t][ix] << endl; }
//			debugfile << "\n\n\n\n\n\n\n\n\n\n\n\n\n"; }
//
//		debugfile.close();
//	}
//	// End of debugging code
//
//
//	return;
//}


///*! This function never actually gets called as of Mar 7 2013*/
//void Trajectory::Initialize_Single_Debris(double *StateVectorPlusUtcIn, double delta_t_in) {
//	// Assuming incoming state vector to be in ECI...though ECEF might be better
//
//	// Currently for single piece of debris
//	int num_debris_pieces = 1;
//
//	if (StateDot_Option != DEBRIS) {
//		cerr << "FATAL ERROR!  Trying to initialize debris, but you chose StateDot = " << StateDot_Option << endl;
//		exit(666); }
//
//	time_start = 0.0;
//	num_per_batch = 1;		// not currently worrying about multiple uncertain runs
//	delta_t = delta_t_in;
//
//	// Unsure how long it will take for debris to hit the ground
//	//  Should not take more than 24 hours!!!
//	int debris_buffer_size = (INTxx) ceil(24*3600/delta_t);
//
//	// Vector of buffer length, where each entry is state vector at that timestep
//	std::vector<double> ZeroVecStateSize(7,0.);
//	State_Vector_Buffer_Vec.assign(debris_buffer_size, ZeroVecStateSize);
//
//	// Load a vector with all the possible times we're willing to consider (spans the buffer)
//	StateTimes = new double [debris_buffer_size];
//	for (int i = 0; i < debris_buffer_size; i++) {
//		StateTimes[i] = time_start + i*delta_t; }
//
//	cout << "Time0 = " << StateTimes[0] << "   Time66 = " << StateTimes[66] << endl;
//
//	//Allocate the storage - NOTE this is NOT a Boost vector, it's std::vector
//	State_Vector_Storage_Vec.assign(num_debris_pieces,std::vector< std::vector<double> >());
//
//
//	//~~~~~~~~~~~~~~~~ This Is the beginning of the propagation! ~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//	//Initial Conditions in ECI
//	ColVec o_R_s__eci(3,0.);
//	ColVec eci_V_s__eci(3,0.);
//
//	// This is just kind of being anal
//	o_R_s__eci(0) = StateVectorPlusUtcIn[0];
//	o_R_s__eci(1) = StateVectorPlusUtcIn[1];
//	o_R_s__eci(2) = StateVectorPlusUtcIn[2];
//	eci_V_s__eci(0) = StateVectorPlusUtcIn[3];
//	eci_V_s__eci(1) = StateVectorPlusUtcIn[4];
//	eci_V_s__eci(2) = StateVectorPlusUtcIn[5];
//	Debris_Mass0 = StateVectorPlusUtcIn[6];
//	Initial_UTC = StateVectorPlusUtcIn[7];		//not actually launch, but breakup time in this case
//
//	// Load up the initial state vector for propagation
//	ColVec state0(state_vec_size);	//[r_vec v_vec mass mass_dot]
//	state0(0) = o_R_s__eci(0);
//	state0(1) = o_R_s__eci(1);
//	state0(2) = o_R_s__eci(2);
//	state0(3) = eci_V_s__eci(0);
//	state0(4) = eci_V_s__eci(1);
//	state0(5) = eci_V_s__eci(2);
//	state0(6) = Debris_Mass0;
//
////	cout << "initial debris vector = " << state0 << endl;
//
//	Initialize_density_profiles(1);		// not really worrying about this at the moment
//
//	timer stopwatch;
//
//	stopwatch.start();
//
//	int piece_number = 0;
////	for (int i = 0; i < num_per_batch; i++) {
//
//	load_density_profile(1);
//
//	int IX = gsl_integrate_debris(state0);
//
//	State_Vector_Storage_Vec[piece_number].assign(IX,ZeroVecStateSize);
//	for (int i = 0; i < IX; i++) {
//		State_Vector_Storage_Vec[piece_number][i] = State_Vector_Buffer_Vec[i]; }
//
//
////	}
//	stopwatch.stop();
//
//
//	//	cout << "that took " << difftime(end_time,start_time) << " seconds" << endl;
//	cout << "in propagate(), that took " << stopwatch.how_long() << " seconds" << endl;
//
//
//
//	write_debris_to_file();
//
//
//
//
//
//
//
//
//
//
//
//
////	// This will have to be allocated on the fly to get the right number of time steps for each piece
////	for (int i = 0; i < num_per_batch; i++) {
////		State_Vector_Storage_Vec[i].assign(num_time_steps,ZeroVecStateSize); }
//
//
//
//
//
//	return;
//}



//// Simulates left/right/hot/cold envelopes given nominal trajectory points
//// Define ellipse that scales with altitude
//// The first point at every time step is the nominal trajectory
//void Trajectory::SimLeftRightHotCold(int fidelity, double semiMajorPercent, double ecc) {
//
//    int num_to_write = 1;
//
//    std::vector<std::vector<std::vector<double> > > RocketLatLonAltStorage = TransformToLatLonAlt(num_to_write);
//
//    //read_single_trajectory_from_file(nominalFileName);
//
//	//Want semi-minor axis along direction of lateral motion, so find direction of motion
//	double LatInitial = RocketLatLonAltStorage[0][0][0];
//	double LonInitial = RocketLatLonAltStorage[0][0][1];
//	double LatFinal = RocketLatLonAltStorage[0][RocketNumTimeSteps[0]-1][0];    //Using the location at stage separation
//	double LonFinal = RocketLatLonAltStorage[0][RocketNumTimeSteps[0]-1][1];
//
//
//	double dy = LatFinal*get_local_earth_radius(LatFinal) - LatInitial*get_local_earth_radius(LatInitial);
//	double dx = LonFinal*get_local_earth_radius(LatFinal)- LonInitial*get_local_earth_radius(LatInitial);
//	double rotAngle = atan2(dy, dx) + PI/2.;	//adding pi/2 to align it with semi-minor
//
//	// higher fidelity means more points used to make the tube
//	double delta_angle = 2.*PI/(fidelity);
//
//
//	// Allocate the memory
//	int num_points_in_trajectory = (INTxx) RocketLatLonAltStorage[0].size();
//	Single_Trajectory_Tube.assign(num_points_in_trajectory, std::vector<Point>());
//	Point EmptyPoint;
//
//
//
//	for (int cur_pt = 0; cur_pt < num_points_in_trajectory; cur_pt++){
//		Single_Trajectory_Tube[cur_pt].assign(fidelity+1,EmptyPoint);
//
//		double gdlat0 = RocketLatLonAltStorage[0][cur_pt][0];
//		double lon0 = RocketLatLonAltStorage[0][cur_pt][1];
//		double alt0 = RocketLatLonAltStorage[0][cur_pt][2];
//		double local_R = get_local_earth_radius(gdlat0);
//
//		// Anticipate some problems here
//		if (alt0 == 0) {
//			alt0 = 0.001; //Set at 1m above ground to avoid division by zero
//		}
//		else if (alt0 < 0) {
//			cout << "ERROR: should not have al0 negative in here, come fix this!!!!~~~~~~~~~" << endl;
//		}
//
//		// Semi-major and Semi-minor axes
//		double a = semiMajorPercent*alt0;
//		double b = a*sqrt(1-ecc*ecc);
//
//		// Location of nominal trajectory on the ground
//		double x0 = R_equator*lon0;
//		double y0 = local_R*gdlat0;
//
//
//
//		//Store the nominal point
//		Single_Trajectory_Tube[cur_pt][0].set_x(x0);
//		Single_Trajectory_Tube[cur_pt][0].set_y(y0);
//		Single_Trajectory_Tube[cur_pt][0].set_z(alt0);
//		Single_Trajectory_Tube[cur_pt][0].set_R_local(local_R);
//
//		//Calculate the envelope points
//		for (int ix = 0; ix < fidelity; ix++) {
//			double r = a*b/sqrt( pow(b*cos(ix*delta_angle),2) + pow(a*sin(ix*delta_angle),2));	//radius from center
//			double xi = r*cos(ix*delta_angle);					// ellipse with semi-major along x-axis
//			double yi = r*sin(ix*delta_angle);
//			double xp = xi*cos(rotAngle) - yi*sin(rotAngle);	// rotate so semi-minor along direction of travel
//			double yp = xi*sin(rotAngle) + yi*cos(rotAngle);
//			double x = x0 + xp;									// offset ellipse around nominal point
//			double y = y0 + yp;
//
//			if (isnan(x)) {
//				cout << "NANS!!!!" << endl;
//			}
//
//			//Store the envelope points
//			Single_Trajectory_Tube[cur_pt][ix+1].set_x(x);
//			Single_Trajectory_Tube[cur_pt][ix+1].set_y(y);
//			Single_Trajectory_Tube[cur_pt][ix+1].set_z(alt0);
//			Single_Trajectory_Tube[cur_pt][ix+1].set_R_local(local_R);
//		} }
//
//	return;
//}


//void Trajectory::write_LRHC_points_instantaneous(char *outFileName, double tstepMin) {
//    // This function writes the LRHC ellipse around the nominal trajectory to a points file so that it may be wrapped up into a footprint
//    // Here, we bin the points according to time
//    //  ASSUMES: Traditional rocket that goes up and never re-enters NAS; will break if you use this for something that re-enters
//    //  ASSUMES: Your first timestep has the rocket within the NAS
//
//	// Do some error checking
//	ofstream outfile;
//	outfile.open(outFileName,ios::out | ios::binary);
//
//	if (outfile.bad()) {
//		cerr << "LRHC FILE FAILED TO OPEN FOR OUTPUT!!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//	// Things we'll need
//	Point temp_pt;
//	std::vector<std::vector<Point> > total_points_at;	//[timestep][points]
//
//    // Figure out which statetimes we need
//    vector <double> StateTimes;
//    if (StateDot_Option == FIRST_STAGE){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == WHOLE_ROCKET){
//        StateTimes = RocketStateTimes;          }
//    else if (StateDot_Option == SUBORBITAL) {
//        StateTimes = SuborbitalStateTimes;      }
//
//	int num_range;
//    // Instead of a solid shape that turns on for a period of time, this will be an ellipse (roughly) propagating upward like a smoke ring
//    // Time binning into minute intervals...may actually make this specific option useless :(
//
//    // We'll want to output this, so save it
//    double startUTC = Initial_UTC;
//
//    //Allocate the vector of points that we're going to write to file
//    int num_points_in_trajectory = (INTxx) Single_Trajectory_Tube.size();   //poorly named, time points, time steps
//
//    double timeInFlight = (StateTimes[num_points_in_trajectory-1] - StateTimes[0])/60.;	//minutes
//    int time_steps_out = (INTxx) ceil(timeInFlight/tstepMin);
//    total_points_at.assign(time_steps_out,std::vector<Point>());
//    num_range = time_steps_out;
//
//    // Collecting the points from a trajectory ellipse / tube
//    if (num_points_in_trajectory > 0) {										//checking to make sure ellipse has been allocated
//        for (int t = 0; t < num_points_in_trajectory; t++) {
//            // NEED SOME LOGIC HERE TO DETERMINE WHAT TIME BIN THE CURRENT RING OF POINTS WILL GO IN
//            int curTimeIX = (INTxx) floor(StateTimes[t]/(60*tstepMin));
//            if (curTimeIX < time_steps_out) {
//                for (int ix = 0; ix < Single_Trajectory_Tube[t].size(); ix++) {
//                    total_points_at[curTimeIX].push_back(Single_Trajectory_Tube[t][ix]); } }
//            else {
//                cerr << "SOMETHING IS WRONG IN TIME-BINNING THE LRHC ELLIPSE!!!~~~~~~~~~~~~~~~~~" << endl;} } }
//    else {
//        cout << "You do not appear to have generated the LRHC boundary~~~~~~~~~~~~~~~~~~~~~~~" << endl; }
//
//
//
//	// Begin writing to file
//	outfile.write((char *) &num_range, sizeof(num_range));
//	outfile.write((char *) &startUTC, sizeof(startUTC));
//	outfile.write((char *) &tstepMin, sizeof(tstepMin));
//
//	int num_points_here;
//	for (int t = 0; t < num_range; t++) {
//		num_points_here = (INTxx) total_points_at[t].size();
//		outfile.write((char *) &num_points_here, sizeof(num_points_here)); }
//
//	for (int t = 0; t < num_range; t++) {
//		outfile.write((char *) &total_points_at[t][0], total_points_at[t].size()*sizeof(Point)); }
//
//	outfile.close();
//
//	// Used for debugging purposes only
//	bool debug = false;
//	if (debug) {
//		string debugFileName("GeneratedFiles/LRHC_PLAINTEXT_DEBUG.txt");
//		ofstream debugfile;
//		debugfile.open(debugFileName.c_str(), ios::out);
//
//		debugfile << num_range << endl << startUTC << endl << tstepMin << endl;
//
//		int num_points_here;
//		for (int t = 0; t < num_range; t++) {
//			num_points_here = (INTxx) total_points_at[t].size();
//			for (int ix = 0; ix < num_points_here; ix++) {
//				debugfile << total_points_at[t][ix] << endl; }
//			debugfile << "\n\n\n\n\n\n\n\n\n\n\n\n\n"; }
//
//		debugfile.close();
//	}
//	// End of debugging code
//
//
//	return;
//}




