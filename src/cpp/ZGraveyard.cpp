

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




/* =============== Old Cumulative FAA , probably not useful anymore because it's been replaced =========================
 
double SkyGrid::generateAllPoints_CumulativeFAA(double thresh, int whichProb){
    //arefMeanList is no longer used.  Keep track of area stats in each cell
    //numberOfPiecesMean is still important because it tells us how many pieces to are in each coefficient group

    //    double MaxIndividualProbability = 0.;  // EV of collision integrated over all timesteps for the given tfail

    // Since this value, for a given converged P_impact distribution, will vary dramatically as bin sizes are changed,
    //      it's really a pretty useless piece of information to return.  Will divide by the correction coefficient
    //      to scale it to be the "maxAllowableProbability per 4nm^2" value
    double maxAllowableProbability = 0.;
//    cout << "CUMULATIVE FAA" << endl;
    //    int numZBins = (INTxx) ceil((NASkm - 0)/zBinHeight);
    double xref = XREF;
    double yref = YREF;
    double zref = ZREF;

    // Nhere
    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh

    // Clear this before we load it up
    all_points_total.clear();
    all_points_total.assign(NumTimestepsHere, vector<Point>());

    generateSpatialProbability(whichProb);
    // This is a map that contains the probabilities of interest (as specified by whichProb) that have been
    //  summed over all time and over all debIX.  Thus, it represents that total probability of, for instance, impact
    //  at a given point in xyz space.
    //    map<int, map<int, map<int,double> > > SpatialProbabilty;

//    // Iterators for the probability grid
//    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
//    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
//    map<int, map<int, map<int,binData> > >::iterator it_x;
//    map<int, map<int, binData> >::iterator it_y;
//    map<int, binData>::iterator it_ID;
//
//    double thisProb = 0;
//    // Have to run through all of the timesteps first to figure out the total probabilities at each timestep
//    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//
//        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
//            int zindex = it_z->first;
//
//            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//
//                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//
//                    it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin();    //Find the first debris index, which is where we store the impact probability
//                    int firstIX = it_ID->first;
//
//                    // Options are PROB_IMPACT, PROB_CASUALTY, PROB_CATASTROPHE
//                    if (whichProb == PROB_IMPACT){
//                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probImpact;}
//                    else if (whichProb == PROB_CASUALTY){
//                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCasualty;}
//                    else if (whichProb == PROB_CATASTROPHE){
//                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCatastrophe;}
//                    else {
//                        cout << "You chose a bad probability option.  Exiting\n";
//                        exit(-99);
//                    }
//
//                    SpatialProbabilty[zindex][xindex][yindex] += thisProb;
//
//                } } }   // Ends loop back to it_z
//    }

    // The thresholds given in the regulations, combined with how the FAA currently computes them, makes them
    //      pointlessly ill-defined because they neither specify what "cumulative" means nor the grid size for the threshold.
    //      So the FAA says: look at grid cell over all time, if the probability of interest exceeds the threshold, then block off that cell.
    //          That's their idea of "cumulative".  I have not seen this definition written anywhere, it has been simply told to me.
    //          This leaves open the possibility of simply shrinking the grid cells until they are so small that none of them exceed
    //          the threshold, then you need not block off any space.  Obviously a bug.
    //          They make some assumptions about aircraft density: specifically one 747 flying in every 4nm^2 square.
    //          I'll use this to update the threshold to be consistent across all cell volumes by making it a threshold DENSITY
    //          This means, also, that using the threshDensity will be MORE CONSERVATIVE than the threshold itself would imply.
    //
    //          We're going to think of this as scaling the threshold
    int numZBins = (INTxx) ceil((NASkm - 0)/zBinHeight);
    double probAircraftPresentInCell = uniformProbabilityValue / numZBins;
    double assumedArea = 13.72;     //4nm^2 in km^2
    double coeff = xBinLength*yBinLength / (assumedArea * numZBins);

    // will eventually be a straight volume comparison: xBinLength*yBinLength*zBinHeight / (assumedArea * NASkm)

    // Now run through one more time to see which points are over the threshold
    map<int, map<int, map<int,double> > >::iterator zit;
    map<int, map<int,double> >::iterator xit;
    map<int,double>::iterator yit;

    vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
    vector<double> temp4Vec;
    temp4Vec.assign(4,0.);

    map<int,int> PointsAtThisLevel;

    for (zit = SpatialProbabilty.begin(); zit != SpatialProbabilty.end(); ++zit){
        int zindex = zit->first;

        for (xit = SpatialProbabilty[zindex].begin(); xit != SpatialProbabilty[zindex].end(); ++xit){
            int xindex = xit->first;

            for (yit = SpatialProbabilty[zindex][xindex].begin(); yit != SpatialProbabilty[zindex][xindex].end(); ++yit){
                int yindex = yit->first;

                // No matter how small the cell area is, assume an aircraft is present and ask "what's the probability of impact?"
                //      Dividing by probAircraftPresentInCell is the "assume an aircraft is present".
                double curProb = SpatialProbabilty[zindex][xindex][yindex] / probAircraftPresentInCell;

                if (curProb >= (thresh * coeff)) {

                    double lowerleftX = xref + xindex*xBinLength;   // Note that xindex is most likely negative
                    double lowerleftY = yref + yindex*yBinLength;   // Note that yindex may be positive or negative
                    double lowerleftZ = zref + zindex*zBinHeight;

                    temp4Vec[0] = lowerleftX;
                    temp4Vec[1] = lowerleftY;
                    temp4Vec[2] = lowerleftZ;
                    temp4Vec[3] = curProb;

                    PointsAtThisLevel[zindex] += 1;
                    ProbabilityHere.push_back(temp4Vec);    }
                else {
                    maxAllowableProbability = std::max(maxAllowableProbability, curProb);
                }

            } } }   // Ends loop back to it_z


    // Save the points
    //      If an xyz point winds up over the limit, then block it off for the whole mission
    //      In a sense, this means there is only one relevent timestep, but let's keep the delta_t unchanged for consistency
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0

        // DEBUGGING
        ProbabilityTotalStorage[tx] = ProbabilityHere;
        stopIXStorage[tx] = -1;       // not relevant here

        // Save the points in the all_points vector
        loadRemainingPointsIntoAllPoints(tx, ProbabilityHere);

    }

//    // Print out the PointsAtThisLevel as a check
//    int numZBins = (INTxx) ceil((NASkm - 0)/zBinHeight);
//    for (int ix = 0; ix < numZBins; ix++){
//        cout << "PointsAtThisLevel[" << ix << "] = " << PointsAtThisLevel[ix] << endl;
//    }


    return (maxAllowableProbability/coeff);
}
=================== END of old FAA cumulative stuff ========================== */





/* =============== Old Cumulative TJC , probably not useful ever =========================

double generateAllPoints_CumulativeTJC(double thresh, int whichProb);

 
// I don't think this is useful anymore.  It's a product of a time when I understood the problem less well.
double SkyGrid::generateAllPoints_CumulativeTJC(double thresh, int whichProb){
    //arefMeanList is no longer used.  Keep track of area stats in each cell
    //numberOfPiecesMean is still important because it tells us how many pieces to are in each coefficient group


//    double MaxIndividualProbability = 0.;  // EV of collision integrated over all timesteps for the given tfail

    double ProbabilityOfImpactSum = 0.;

    //    int numZBins = (INTxx) ceil((NASkm - 0)/zBinHeight);
    double xref = XREF;
    double yref = YREF;
    double zref = ZREF;

    // Nhere
    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh

    // Clear this before we load it up
    all_points_total.clear();
    all_points_total.assign(NumTimestepsHere, vector<Point>());

    // Iterators for the probability grid
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;

    double thisProb = 1.;

    int tfail = 0;
    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {

        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
        if (it_time==ProbabilityMapDebIX.begin()){ tfail = tx; }

        vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
        vector<double> temp4Vec;
        temp4Vec.assign(4,0.);

        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
            int zindex = it_z->first;
            //cout << "  " << xindex << endl;

            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                int xindex = it_x->first;
                //cout << "    " << yindex << endl;

                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                    int yindex = it_y->first;

                    it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin();    //Find the first debris index, which is where we store the impact probability
                    int firstIX = it_ID->first;

                    // Options are PROB_IMPACT, PROB_CASUALTY, PROB_CATASTROPHE
                    if (whichProb == PROB_IMPACT){
                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probImpact;}
                    else if (whichProb == PROB_CASUALTY){
                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCasualty;}
                    else if (whichProb == PROB_CATASTROPHE){
                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCatastrophe;}
                    else {
                        cout << "You chose a bad probability option.  Exiting\n";
                        exit(-99);
                    }


                    double lowerleftX = xref + xindex*xBinLength;   // Note that xindex is most likely negative
                    double lowerleftY = yref + yindex*yBinLength;   // Note that yindex may be positive or negative
                    double lowerleftZ = zref + zindex*zBinHeight;

                    temp4Vec[0] = lowerleftX;
                    temp4Vec[1] = lowerleftY;
                    temp4Vec[2] = lowerleftZ;
                    temp4Vec[3] = thisProb;

                    ProbabilityHere.push_back(temp4Vec);

                    if (temp4Vec[3] < 0){
                        // THIS ONE IS EXTRA IMPORTANT!!!
                        // http://stackoverflow.com/questions/8480640/how-to-throw-a-c-exception
                        // cout << "[" << tx << ", " << xindex << ", " << yindex << ", " << zindex << "] = " << temp4Vec[3] << "  YOU WENT NEGATIVE!!!\n\n";
                    }
                } } }   // Ends loop back to it_z

        // Since the grid cell areas are the same AND assuming uniform aircraft density, sorting the CDF gives same ranking as sorting the EsubC values
        // using myfunction: compares the third index (the probability) of each vector passed to it
        // these should be sorted from highest to lowest probability
        std::sort(ProbabilityHere.begin(), ProbabilityHere.end(), compareXYZ);

        // ~~~~~~~~ TJC Conception of Cumulative Probability Protection ~~~~~~~~~~~~~
        int stopIX = ProbabilityHere.size();      // index in reverse (lowest prob first)
        double probSum = 0.;            // EsubC = probSum * ProbOfAirplane * ExpectedValueNumDebrisHere
        double probLevel = thresh / (NumTimestepsHere - tfail + 1);     // SHOULD SUBTRACT DENOMINATOR BY TFAIL (first tx)!!!!

        // Figure out where the last index is
        while (probSum < probLevel){
            stopIX--;                   // Subtract 1 before accessing vector for first time
            probSum += ProbabilityHere[stopIX][3];
            if (stopIX == 0){ break; }
        }

        // Two options here that make stopIX == 0:
        //  1.) because there are actually no points that need to be excluded
        //  2.) the very last point (perhaps there is only one point total) puts you over the threshold

        // If that last index put us over the limit, add one back.
        if (probSum > probLevel){
            probSum -= ProbabilityHere[stopIX][3];
            stopIX++;
        }

        ProbabilityOfImpactSum += probSum;

        // DEBUGGING
        ProbabilityTotalStorage[tx] = ProbabilityHere;
        stopIXStorage[tx] = stopIX;

        // resize the vector and return it
        ProbabilityHere.resize(stopIX);

        // Save the points in the all_points vector
        loadRemainingPointsIntoAllPoints(tx, ProbabilityHere);
    }

    return ProbabilityOfImpactSum;
}
=================== END of TJC cumulative stuff ========================== */


















