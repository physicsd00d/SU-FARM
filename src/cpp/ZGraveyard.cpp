

/* ================== Stuff for creating aircraft density maps =========================
vector<double> createEmptyAircraftDensityMap();
void populateAircraftDensityMap(void *densityMapArray, int numElements);
bool usingMITDensityMap;
double uniformProbabilityValue;


// ACDensityMap doesn't appear to get used anywhere and it's not actually a good idea for mission planning.
// TODO: REMOVE EVERYTHING that has to do with this.
vector<double> SkyGrid::createEmptyAircraftDensityMap(){

//    map<int, map<int, map<int, double> > > ACDensityMap;

    {   // Scoping this so the iterators fall out of scope immediately afterwards
        map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
        map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
        map<int, map<int, map<int,binData> > >::iterator it_x;
        map<int, map<int, binData> >::iterator it_y;

        // First loop through the indices and find the unique xy cells (omits t and z)
        for (it_time=ProbabilityMapDebIX.begin(); it_time!=ProbabilityMapDebIX.end(); it_time++) {
            int tx = it_time->first;

            for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
                int zindex = it_z->first;

                for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                    int xindex = it_x->first;

                    for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                        int yindex = it_y->first;

                        ACDensityMap[0][xindex][yindex] = -1.;
                    } } } }
    }

    // Now load up the vector of those cells
    vector<double> flatLatLons;

    map<int, map<int, double> >::iterator it_x;
    map<int, double>::iterator it_y;

    Point NW;
    Point SE;

    for (it_x = ACDensityMap[0].begin(); it_x != ACDensityMap[0].end(); ++it_x){
        int xindex = it_x->first;

        double leftX    = XREF + xindex*xBinLength;
        double rightX   = XREF + (xindex+1)*xBinLength;
        for (it_y = ACDensityMap[0][xindex].begin(); it_y != ACDensityMap[0][xindex].end(); ++it_y){
            int yindex = it_y->first;

            double bottomY  = YREF + yindex*yBinLength;
            double topY     = YREF + (yindex+1)*yBinLength;

            NW.set_xyz(leftX, topY, 0);
            SE.set_xyz(rightX, bottomY, 0);

            flatLatLons.push_back(NW.get_gdLatDeg());
            flatLatLons.push_back(NW.get_lonDeg());
            flatLatLons.push_back(SE.get_gdLatDeg());
            flatLatLons.push_back(SE.get_lonDeg());
        }
    }

//    // Convert vector to double array because i think python can handle that
//    int numLatLonPts = flatLatLons.size();
//    double *flatLatLonArray = new double [numLatLonPts];
//    memcpy((void *)flatLatLonArray, (void *)&flatLatLons[0], numLatLonPts*sizeof(double));
//
//    for (int ix = 0; ix < numLatLonPts; ix++){
//        cout << "ix = " << ix << ",  pt = " << flatLatLonArray[ix] << endl;
//    }


    return flatLatLons;
}


// ACDensityMap doesn't appear to get used anywhere and it's not actually a good idea for mission planning.
// TODO: REMOVE EVERYTHING that has to do with this.
void SkyGrid::populateAircraftDensityMap(void *densityMapArray, int numElements){

    if (numElements > 0){
        // Make a record that we're using a custom map
        usingMITDensityMap = true;

        // leave uniformProbabilityValue unset at -1

        vector<double> probArray ((double*)densityMapArray, (double*)densityMapArray + numElements);

        map<int, map<int, double> >::iterator it_x;
        map<int, double>::iterator it_y;

        int counter = 0;
        for (it_x = ACDensityMap[0].begin(); it_x != ACDensityMap[0].end(); ++it_x){
            int xindex = it_x->first;

            for (it_y = ACDensityMap[0][xindex].begin(); it_y != ACDensityMap[0][xindex].end(); ++it_y){
                int yindex = it_y->first;

                ACDensityMap[0][xindex][yindex] = probArray[counter];
                counter++; } }
    } else if (numElements == -1){
        // Make a record that we're using a uniform map
        usingMITDensityMap = false;

        // Populate the ACDensityMap with a constant value specified by the user
        double probValue = ((double*)densityMapArray)[0];
        uniformProbabilityValue = probValue;

        map<int, map<int, double> >::iterator it_x;
        map<int, double>::iterator it_y;

        int counter = 0;
        for (it_x = ACDensityMap[0].begin(); it_x != ACDensityMap[0].end(); ++it_x){
            int xindex = it_x->first;

            for (it_y = ACDensityMap[0][xindex].begin(); it_y != ACDensityMap[0][xindex].end(); ++it_y){
                int yindex = it_y->first;

                ACDensityMap[0][xindex][yindex] = probValue;
                counter++; } }

    }

    return;
}
====================== End of aircraft density map stuff =================== */

























