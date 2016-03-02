//
//  SkyGrid.cpp
//  Prop3Dof
//
//  Created by Thomas Colvin on 9/25/13.
//  Copyright (c) 2013 Thomas Colvin. All rights reserved.
//

#include "SkyGrid.h"

// Constructor that takes a points file
SkyGrid::SkyGrid(string CapeLrhcFile, double xBinLength_in, double yBinLength_in, double zBinHeight_in)
        : PointCloud(CapeLrhcFile)
{
    xBinLength = xBinLength_in;
    yBinLength = yBinLength_in;
    zBinHeight = zBinHeight_in;
    
//    halfBufferCells = ceil(bufferDist/(std::min(xBinLength,yBinLength)));
    
//    GridTheSky(all_points_total);
    GridTheSky();
    
    int numTimeSteps = all_points_total.size();
    
    isProbability = false;
    doneASH = false;
    hazardProbabilitiesGenerated = false;
//    uniformProbabilityValue = -1;

    
    return;
}

// Constructor that takes a pointer to a PointCloud object
SkyGrid::SkyGrid(PointCloud *newCloud, double xBinLength_in, double yBinLength_in, double zBinHeight_in)
        : PointCloud(newCloud)
{
    xBinLength = xBinLength_in;
    yBinLength = yBinLength_in;
    zBinHeight = zBinHeight_in;
    
//    // Because the prob density of being in a cell depends on the cell volume (thus zBinHeight)
//    //      and because we throw away points that are above the NAS...we need to cap the zBinHeight
//    //      at the top of the NAS otherwise we've got a bug that would allow you to shrink your probability
//    //      of finding debris in a cell by simply letting zBinHeight go to infinity.
//    if (zBinHeight > NASkm){ zBinHeight = NASkm; }
    
//    newCloud->PrintAllPoints();
    
    GridTheSky();
        
    int numTimeSteps = all_points_total.size();
    
    isProbability = false;
    doneASH = false;
    hazardProbabilitiesGenerated = false;
//    uniformProbabilityValue = -1;


//    cout << "addy of newCloud = " << newCloud << "   " << &newCloud << endl;
    
    return;
}
 
  
// Functions for probabilistically handling footprints 
void SkyGrid::GridTheSky(){
    //total_points_at[tstep][piece of debris] 
    
    // Pick a reference point (x,y,z) = (0,0,0) and i suppose all_points_UTC for now
    double xref = XREF;
    double yref = YREF;
    double zref = ZREF;
    
//    map<int, int> countDebrisIX;

    
//    vector<vector<Point> > total_points_at = getAllPoints();
    
//    int numGoodPts = 0;
    
    // Bin the points
    unsigned long numTimeSteps = all_points_total.size();
//    cout << "numTimeSteps = " << numTimeSteps << endl;

    for (int tx = 0; tx < numTimeSteps; tx++){

        for (int px = 0; px < getNumPointsAtTstep(tx); px++){
//            Point curPoint = total_points_at[tx][px];
            Point curPoint = all_points_total[tx][px];
            
            int curID = curPoint.get_id();
            int xindex = floor((curPoint.get_x() - xref)/xBinLength);
            int yindex = floor((curPoint.get_y() - yref)/yBinLength);
            int zindex = floor((curPoint.get_z() - zref)/zBinHeight);
            double thisMass = curPoint.mass_kg;
            double thisArea = curPoint.area_km2;
            
//            countDebrisIX[curID] += 1;  // Mainly i just want to know how many curIDs there are, but might as well keep track of total num pts

            double Vx = curPoint.Vx_km_s;
            double Vy = curPoint.Vy_km_s;
            double Vz = curPoint.Vz_km_s;
            
            // Adjust zindex
            // If z is negative, place it on the ground, if z is too high, ignore it
            bool zgood = true;
            if (curPoint.get_z() > NASkm) {
                zindex = -2;
                xindex = 0;
                yindex = 0;
                Vx = 0.;
                Vy = 0.;
                Vy = 0.;
                thisMass = 0.;
                thisArea = 0.;
            } else if (curPoint.get_z() < 0) {
                zindex = -1;
                xindex = 0;
                yindex = 0;
                Vx = 0.;
                Vy = 0.;
                Vy = 0.;
                thisMass = 0.;
                thisArea = 0.;
            }

//            // If z is ever out of bounds, either too high or already landed, save it anyways.
//            bool zgood = true; // Will always be good now.
//            if ((curPoint.get_z() > NASkm) || (curPoint.get_z() < 0)) {
//                zindex = -1;
//                xindex = 0;
//                yindex = 0;
//                Vx = 0.;
//                Vy = 0.;
//                Vy = 0.;
//                thisMass = 0.;
//                thisArea = 0.;
//            }
            
//            printf("[%d][%d][%d]\n", tx, curID, zindex);
            
            // Might also get passed in spurious points (too high to care about), skip them
            //   only add the points that are currently within the NAS
            if (zgood) {
                // Keep track of this because it will be needed for checking that the probabilities are consistent
                numPtsLessThanReactionInNAS[curID] += 1;
                
                // Compute the magnitude of the velocity
                double velNorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
                
                // Get the binData for the current cell
                binData curData = GridMapDebIX[tx][curID][zindex][xindex][yindex];
                double countsHere = curData.probDebris;
                
                // Update the velocity info
                curData.avgVel = ((countsHere*(curData.avgVel) + velNorm)/(countsHere + 1));
                curData.maxVel = std::max(velNorm, curData.maxVel);
                curData.minVel = std::min(velNorm, curData.minVel);
                
                // Update the mass info
                curData.avgMass = ((countsHere*(curData.avgMass) + thisMass)/(countsHere + 1));
                curData.maxMass = std::max(thisMass, curData.maxMass);
                curData.minMass = std::min(thisMass, curData.minMass);
                
                // Update the area info
                curData.avgArea = ((countsHere*(curData.avgArea) + thisArea)/(countsHere + 1));
                curData.maxArea = std::max(thisArea, curData.maxArea);
                curData.minArea = std::min(thisArea, curData.minArea);
                
                // Update the count
                curData.probDebris += 1;
                
                totalNumPtsAtStepMapDebIX[tx][curID] += 1;      // This is not even all of the points that are below the reaction time
                                                                //  I don't think this should get used for anything, it's not meaningful
                                                                //  given the statistics I'm trying to gather.  This would only benefit me
                                                                //  if i wanted to find the probability of debris GIVEN (under reaction time
                                                                //  AND inside the NAS).
                                                                //  DON'T USE THIS ANYWHERE!  THAT'S AN ORDER!
                                                                //  CALM DOWN!  This might be good as a probability check
                
                
                
                
                GridMapDebIX[tx][curID][zindex][xindex][yindex] = curData;
            }
        }
        
//        cout << "t = " << tx << ",  maxVel = " << maxVel << ",  maxAvgVel = " << maxAvgVel <<  endl;
    }
     
    
//    // Dump some stats
//    map<int, map<int, int> >::iterator it;
//    map<int, int>::iterator id;
//    cout << "what?" << endl;
//    for (it=totalNumPtsAtStepMapDebIX.begin(); it!=totalNumPtsAtStepMapDebIX.end(); it++) {
//        printf("t = %d\n", it->first);
//        for (id = it->second.begin(); id != it->second.end(); ++id){
//            printf("     debIX = %d,   num = %d\n", id->first, id->second);
//        } }
    
//    cout << "numGoodPts = " << numGoodPts << endl;
    
    return;
}




void SkyGrid::PythonDebrisIntoGrid(PointCloud *in){
    
    // Check that you can still use it
    if (isProbability){
        cout << "ERROR, youre trying to use PythonDebrisIntoGrid after converting to probabilities\n\n";
        exit(-13);
    }
    
//    double inTstepMinutes  = in->getDeltaT();
    double inTstepSeconds  = in->getDeltaT();
    double inInitialUTC    = in->getInitialUTC();
//    int inNumRange         = in->getPointsRange();
    
//    double curTstepMinutes  = getDeltaT();
    double curTstepSeconds  = getDeltaT();
    double curInitialUTC    = getInitialUTC();
//    int curNumRange         = Grid.size();
    
    // Checks that the timestep sizes are the same
    double deltaTsec = 0;
    if (inTstepSeconds != curTstepSeconds){
        cout << "ERROR!!  TIMESTEPS DON'T MATCH!!!  RETURNING FROM FUNCTION WITHOUT DOING ANYTHING!!!!!";
        exit(-10);
    } else {
        deltaTsec = curTstepSeconds;
    }
    
//    // Check the time range of the new points, may need to extend the time vector in either direction
//	double curFinalUTC = curInitialUTC + curTstepMinutes*curNumRange/(60*24);
//	double inFinalUTC = inInitialUTC + inTstepMinutes * inNumRange/(60.*24.);
    
    // Determines if the incoming points overstep the time bounds in either direction
    bool isLeading = (inInitialUTC < curInitialUTC);
    
//    // As long as this is not zero, the grid will be recalculated.
//    int leadingTimeSteps = 0;
    
    if (isLeading){
        cout << "For debugging reasons, i'm going to exit the program if there are any leading timesteps...sorry....\n";
        exit(-12);
    }
    
    // NEED TO WORRY ABOUT IF THE INCOMING POINTS START ***AFTER*** THE CURRENT POINTS
    int incomingOffset = max((int) round((inInitialUTC - curInitialUTC)*24*3600./deltaTsec), 0);
    if (incomingOffset != 0){
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         icomingOffset = " << incomingOffset << endl;
        exit(-9);
    }

    // Copy the allPoints from the incoming cloud to the current cloud
    //    FOR SAFETY: Should probably make sure that I've already copied over the other incoming values that I would want to keep
    all_points_total = in->getAllPoints();
    
    GridTheSky();
    
    return;
}


// This is used for ordering the probabilities in the histograms
bool compareXYZ (vector<double> i, vector<double> j) { return (i[3] > j[3]); }









//double SkyGrid::generateAllPointsFromASH(vector<int> pieces, vector<double> areas, int numDebIX, double thresh, double pFail){
//    
//    if (!doneASH){
//        // You didn't do the ASH yet!
//        cout << "ERROR!!!!  You didn't do the ASH yet!" << endl;
//        return(-100);
//        
//        // need to learn how to raise errors
//    }
//    
//    // Options are PROB_IMPACT, PROB_CASUALTY, PROB_CATASTROPHE
//    int whichProb = PROB_IMPACT;
//    generateProbabilityOfImpact(pieces, pFail);
//    
//    
////    double EV_strike = generateAllPoints_CumulativeTJC(thresh, whichProb);
//    
//    double EV_strike = generateAllPoints_CumulativeFAA(thresh, whichProb);
//    
////    PrintAllPoints();
//    
//    // I don't think i included the whichProb into this one, so watch out
//    // double EV_strike = generateAllPointsFromProbability(pieces, areas, thresh, pFail);
//    
//    return EV_strike;
//}





/*! Despite its name, it actually generates all three possible probabilities:
 * * Prob Impact
 * * Prob Casualty
 * * Prob Catastrophe
 *
 * Note that it partially destroys the ProbabilityMapDebIX structure by writing the values
 *  of the probabilities to the first index of the debris category...so don't try to use this again.
 *  This sucks, I know, but it saves some memory and i have other things to do right now.
 *
 * Removing pFail, this calculates \prod_d^D (1 - \phi_{ij}  \Xi^d_{i \mid jf})
 */
void SkyGrid::generateHazardProbabilities(vector<int> numberOfPiecesMean){
    
    // Check to make sure that this hasn't already been done
    if (not hazardProbabilitiesGenerated){
        // Proceed, but flip the flag so we know we've been here
        hazardProbabilitiesGenerated = true; }
    else{
        cout << "ERROR!!!!  You cannot call generateHazardProbabilities twice!" << endl;
        exit(-100); }
    
    // ============== Estimation of 787 cruising through airspace, this is an OVERESTIMATE!!!  See RCC ===========
    // For a 787
    // typical cruising speed is 900km/h
    double speed787 = 954./3600.;     //km/s      (954km/h = 15.9km/min)
    
    //    double fourNM2 = 13.72;         // 4 (n.m.)^2 * (1.852 km/nm)^2 = 13.72 km^2
    //    double aircraftDensity = 1./fourNM2;       // [prob/km^2] Paul Wilde's assumed aircraft density (1 every 4nm^2)
    
    // (1/13.72) * 5*5/4
    
    //    double cellArea = xBinLength*yBinLength;
    //    double probOfAirplaneInCell = aircraftDensity * cellArea / numZBins;    // Probability of having an aircraft in a cell (divide by number of zbins)
    double ft_2_km = 0.0003048;
    
    // Quasi-following Wilde's Modeling of Risk to Aircraft from Space Vehicle Debris
    //    double span787 = 0.0601401;     //km        (span = 197.31ft  =  0.0601401km)
    //    double length787 = 0.0565587;   //km        (length = 185.56ft = 0.0565587km)
    //    double height787 = 0.010;       //km        ASSUMING a height of 10m, look up the real number later
    //    double d_Airplane_front = sqrt(span787 * height787);
    //    double d_Airplane_top = sqrt(span787 * length787);
    
    // Looking at RCC 321, they say maximum area is 12,000ft for B747.  Based on that, I estimate these
    double d_Airplane_top = sqrt(1018.5 /* m^2 */) * (1e-3);         // This is plan area + length*cabinwidth   [km]
    double d_Airplane_front = sqrt(96.3405 /* m^2 */) * (1e-3);      // This is 12,000ft^2 - TopArea            [km]
    
    //    double volumeToAvoid = span787*(speed787 * delta_t + length787) * height787;
    double cellVolume = xBinLength*yBinLength*zBinHeight;
    double delta_t = getDeltaT();
    
    // Do a little bit of error checking
    double maxProbSingleStrike = 1.;
    double minProbSingleStrike = 0.;
    
    double maxSweptVolume = cellVolume;
    
    // ===========  Start to load up the probabilities at the given time index tx
    //    vector<double> temp4Vec;
    //    temp4Vec.assign(4,0.);
    
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    
    // Find the probability vector at this timestep for all points in the grid.
    //    vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
    
    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
        
        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
        
        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
            int zindex = it_z->first;
            //cout << "  " << zindex << endl;
            
            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                int xindex = it_x->first;
                //cout << "    " << yindex << endl;
                
                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                    int yindex = it_y->first;
                    
                    // ========== This is roughly how Wilde computes the probabilities =============
                    double probNoStrike         = 1.;
                    double probNoCasualty       = 1.;
                    double probNoCatastrophe    = 1.;

                    // Define the index where we'll store info and then check that it's not already being used.
                    if (ProbabilityMapDebIX[tx][zindex][xindex][yindex].count(STORE_IX) > 0){
                        cout << "ERROR: Some debris class has stolen the ID used for recording hazard probabilities" << endl;
                        exit(-2);
                    }
                    
                    for (it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin(); it_ID != ProbabilityMapDebIX[tx][zindex][xindex][yindex].end(); ++it_ID){
                        int curID = it_ID->first;
                        binData PD = (it_ID->second);
                        
                        // Assuming that if we know there is debris in this cell, the location of that debris is uniformly likely to be anywhere in the volume.
                        // Keep in mind that as the grid gets finer, probDebrisHere can get bigger than 1 or arbitrarily large because we're dividing by such a
                        // small number.
                        
                        double probDebrisInCell     = PD.probDebris;
                        double probDensityOfDebris  = 1. / cellVolume;    //Swept volume of debris / cellVolume
                        
                        //                        double probDebrisHere = ((it_ID->second).probDebris)/cellVolume;        // Prob Density
                        
                        //                        double AF = pow(d_Airplane_front + sqrt(PD.avgArea),2);
                        //                        double AT = pow(d_Airplane_top + sqrt(PD.avgArea),2);
                        //                        double probOfSingleStrike = probDebrisHere * (AF*speed787 + AT*PD.avgVel) * delta_t;
                        
                        // Equation takes mass in grams, outputs ft^2, so need to convert to km^2
                        double theta = atan2(speed787, PD.avgVel);  // Sure hope velocity is in km/s
                        double A_Proj = pow(d_Airplane_front, 2)*sin(theta) + pow(d_Airplane_top, 2)*cos(theta);
                        double V_impact = sqrt(pow(speed787, 2) + pow(PD.avgVel,2));
                        
                        double A_Casualty   = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);   // This is the case for pieces over 300g
                        double A_Catastrope = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);
                        
                        if (PD.avgMass < 0.001){
                            // Throw out pieces that are less than 1g.  They pose no danger
                            A_Casualty      = 0.;
                            A_Catastrope    = 0.;
                            A_Proj          = 0.;
                        }
                        else if (PD.avgMass < 0.300){
                            // If the mass is below 300g, these areas are modeled as a function of the piece's mass
                            A_Casualty      = ((0.0085*pow(PD.avgMass * 1e3,2) + 8.5*(PD.avgMass * 1e3) + 200)) * pow(ft_2_km,2);
                            A_Catastrope    = (0.025 * pow(PD.avgMass * 1e3,2) + 4*(PD.avgMass * 1e3)) * pow(ft_2_km,2);
                        }
                        
                        // Prob density integrated over swept volume of airplane
                        //                        double probOfSingleStrike   = probDebrisInCell * probDensityOfDebris * A_Proj * V_impact * delta_t;         // Pure Probability
                        double probOfSingleStrike   = probDebrisInCell * probDensityOfDebris * A_Proj * V_impact * delta_t;         // Pure Probability
                        double probOfCasualty       = probDebrisInCell * probDensityOfDebris * A_Casualty * V_impact * delta_t;
                        double probOfCatastrophe    = probDebrisInCell * probDensityOfDebris * A_Catastrope * V_impact * delta_t;
                        
                        if (probOfSingleStrike > maxProbSingleStrike){
                            printf("ERROR: probOfSingleStrike %E > %E maxProbSingleStrike\n",probOfSingleStrike, maxProbSingleStrike);
                            exit(-90);
                        } else if (probOfSingleStrike < minProbSingleStrike){
                            printf("ERROR: probOfSingleStrike %E > %E minProbSingleStrike\n",probOfSingleStrike, minProbSingleStrike);
                            exit(-90);
                        } else if (A_Proj * V_impact * delta_t > maxSweptVolume){
                            printf("ERROR: Vswept %E > %E maxSweptVolume\n", A_Proj * V_impact * delta_t, maxSweptVolume);
                            exit(-90);
                        }
                        
                        // How different are the projected areas anyways?
                        // printf("Proj %E  --  Cas %E  --  Cat %E\n", A_Proj, A_Casualty, A_Catastrope);
                        
                        double expectedNumPiecesHere        = numberOfPiecesMean[curID];      // This is the number of debris from this debIX in the catalog
                        double probOfNoStrikeFromCurID      = pow(1. - probOfSingleStrike, expectedNumPiecesHere);
                        double probOfNoCasualtyFromCurID    = pow(1. - probOfCasualty, expectedNumPiecesHere);
                        double probOfNoCatastropheFromCurID = pow(1. - probOfCatastrophe, expectedNumPiecesHere);
                        
                        probNoStrike        *= probOfNoStrikeFromCurID;
                        probNoCasualty      *= probOfNoCasualtyFromCurID;
                        probNoCatastrophe   *= probOfNoCatastropheFromCurID;

//                        if ((zindex < 0)){
//                            printf("zindex = %d, probDebrisInCell = %E, Aproj = %E, A_Casualty = %E \n", zindex, probDebrisInCell, A_Proj, A_Casualty);
//                        }
                    }
                    
                    // Error check on LandedIX and TooHighIX.  Should be exactly zero probability of strike, probNoStrike = 1.
                    // This is because both mass and velocity have been set to zero when gridding, so the projected areas will
                    // be zero.  Thus probOfNoStrikeFromCurID = (1.-0.)^numPieces should be basically exactly 1.  This is kind
                    // of a fragile condition, checking for equality with a double, but it really should be exactly 1.0.  If
                    // it's not, then something is wrong.  Maybe should just set it to 1.0 for safety?
                    if ((zindex < 0) && (probNoStrike != 1.)){
                        printf("ERROR generateHazardProbabilities: zindex = %d  vs probNoStrike = %E \n", zindex, probNoStrike);
                        exit(-12);
                    }
                    
                    // Save the probability for this grid cell.  Only put it in the leading debIX.
                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][STORE_IX].probNoImpact         = probNoStrike;
                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][STORE_IX].probNoCasualty       = probNoCasualty;
                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][STORE_IX].probNoCatastrophe    = probNoCatastrophe;
                    
                } } } }
    
    return;
}



// Commented out Feb 2016
///*! Despite its name, it actually generates all three possible probabilities:
// * * Prob Impact
// * * Prob Casualty
// * * Prob Catastrophe
// *
// * Note that it partially destroys the ProbabilityMapDebIX structure by writing the values
// *  of the probabilities to the first index of the debris category...so don't try to use this again.
// *  This sucks, I know, but it saves some memory and i have other things to do right now.
// */
////void SkyGrid::generateProbabilityOfImpact(vector<int> numberOfPiecesMean, double pFail){
//void SkyGrid::generateHazardProbabilities(vector<int> numberOfPiecesMean, double pFail){
//    
//    // Check to make sure that this hasn't already been done
//    if (not hazardProbabilitiesGenerated){
//        // Proceed, but flip the flag so we know we've been here
//        hazardProbabilitiesGenerated = true; }
//    else{
//        cout << "ERROR!!!!  You cannot call generateHazardProbabilities twice!" << endl;
//        exit(-100); }
//    
//    // ============== Estimation of 787 cruising through airspace, this is an OVERESTIMATE!!!  See RCC ===========
//    // For a 787
//    // typical cruising speed is 900km/h
//    double speed787 = 954./3600.;     //km/s      (954km/h = 15.9km/min)
//    
//    //    double fourNM2 = 13.72;         // 4 (n.m.)^2 * (1.852 km/nm)^2 = 13.72 km^2
//    //    double aircraftDensity = 1./fourNM2;       // [prob/km^2] Paul Wilde's assumed aircraft density (1 every 4nm^2)
//    
//    // (1/13.72) * 5*5/4
//    
//    //    double cellArea = xBinLength*yBinLength;
//    //    double probOfAirplaneInCell = aircraftDensity * cellArea / numZBins;    // Probability of having an aircraft in a cell (divide by number of zbins)
//    double ft_2_km = 0.0003048;
//    
//    // Quasi-following Wilde's Modeling of Risk to Aircraft from Space Vehicle Debris
//    //    double span787 = 0.0601401;     //km        (span = 197.31ft  =  0.0601401km)
//    //    double length787 = 0.0565587;   //km        (length = 185.56ft = 0.0565587km)
//    //    double height787 = 0.010;       //km        ASSUMING a height of 10m, look up the real number later
//    //    double d_Airplane_front = sqrt(span787 * height787);
//    //    double d_Airplane_top = sqrt(span787 * length787);
//    
//    // Looking at RCC 321, they say maximum area is 12,000ft for B747.  Based on that, I estimate these
//    double d_Airplane_top = sqrt(1018.5 /* m^2 */) * (1e-3);         // This is plan area + length*cabinwidth   [km]
//    double d_Airplane_front = sqrt(96.3405 /* m^2 */) * (1e-3);      // This is 12,000ft^2 - TopArea            [km]
//    
//    //    double volumeToAvoid = span787*(speed787 * delta_t + length787) * height787;
//    double cellVolume = xBinLength*yBinLength*zBinHeight;
//    double delta_t = getDeltaT();
//    
//    // Do a little bit of error checking
//    double maxProbSingleStrike = 0.;
//    double minProbSingleStrike = 1.;
//    
//    double maxSweptVolume = 0.;
//    
//    // ===========  Start to load up the probabilities at the given time index tx
//    //    vector<double> temp4Vec;
//    //    temp4Vec.assign(4,0.);
//    
//    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
//    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
//    map<int, map<int, map<int,binData> > >::iterator it_x;
//    map<int, map<int, binData> >::iterator it_y;
//    map<int, binData>::iterator it_ID;
//    
//    
//    // Find the probability vector at this timestep for all points in the grid.
//    //    vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
//    
//    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
//        
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//        
//        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
//            int zindex = it_z->first;
//            //cout << "  " << xindex << endl;
//            
//            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//                //cout << "    " << yindex << endl;
//                
//                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//                    
//                    // ========== This is roughly how Wilde computes the probabilities =============
//                    double probNoStrike         = 1.;
//                    double probNoCasualty       = 1.;
//                    double probNoCatastrophe    = 1.;
//                    
//                    int firstIX = -666;
//                    for (it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin(); it_ID != ProbabilityMapDebIX[tx][zindex][xindex][yindex].end(); ++it_ID){
//                        int curID = it_ID->first;
//                        
//                        if (firstIX == -666) { firstIX = curID; }
//                        
//                        binData PD = (it_ID->second);
//                        // Assuming that if we know there is debris in this cell, the location of that debris is uniformly likely to be anywhere in the volume.
//                        // Keep in mind that as the grid gets finer, probDebrisHere can get bigger than 1 or arbitrarily large because we're dividing by such a
//                        // small number.
//                        
//                        double probDebrisInCell     = PD.probDebris;
//                        double probDensityOfDebris  = 1. / cellVolume;    //Swept volume of debris / cellVolume
//                        
////                        double probDebrisHere = ((it_ID->second).probDebris)/cellVolume;        // Prob Density
//                        
//                        //                        double AF = pow(d_Airplane_front + sqrt(PD.avgArea),2);
//                        //                        double AT = pow(d_Airplane_top + sqrt(PD.avgArea),2);
//                        //                        double probOfSingleStrike = probDebrisHere * (AF*speed787 + AT*PD.avgVel) * delta_t;
//                        
//                        // Equation takes mass in grams, outputs ft^2, so need to convert to km^2
//                        double theta = atan2(speed787, PD.avgVel);  // Sure hope velocity is in km/s
//                        double A_Proj = pow(d_Airplane_front, 2)*sin(theta) + pow(d_Airplane_top, 2)*cos(theta);
//                        double V_impact = sqrt(pow(speed787, 2) + pow(PD.avgVel,2));
//                        
//                        double A_Casualty   = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);   // This is the case for pieces over 300g
//                        double A_Catastrope = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);
//                        
//                        if (PD.avgMass < 0.001){
//                            // Throw out pieces that are less than 1g.  They pose no danger
//                            A_Casualty      = 0.;
//                            A_Catastrope    = 0.;
//                        }
//                        else if (PD.avgMass < 0.300){
//                            // If the mass is below 300g, these areas are modeled as a function of the piece's mass
//                            A_Casualty      = ((0.0085*pow(PD.avgMass * 1e3,2) + 8.5*(PD.avgMass * 1e3) + 200)) * pow(ft_2_km,2);
//                            A_Catastrope    = (0.025 * pow(PD.avgMass * 1e3,2) + 4*(PD.avgMass * 1e3)) * pow(ft_2_km,2);
//                        }
//                        
//                        // Prob density integrated over swept volume of airplane
////                        double probOfSingleStrike   = probDebrisInCell * probDensityOfDebris * A_Proj * V_impact * delta_t;         // Pure Probability
//                        double probOfSingleStrike   = probDebrisInCell * probDensityOfDebris * A_Catastrope * V_impact * delta_t;         // Pure Probability
//                        double probOfCasualty       = probDebrisInCell * probDensityOfDebris * A_Casualty * V_impact * delta_t;
//                        double probOfCatastrophe    = probDebrisInCell * probDensityOfDebris * A_Catastrope * V_impact * delta_t;
//                        
//                        maxProbSingleStrike = std::max(maxProbSingleStrike  , probOfSingleStrike); 
//                        minProbSingleStrike = std::min(minProbSingleStrike  , probOfSingleStrike);
//                        maxSweptVolume      = std::max(maxSweptVolume       , A_Proj * V_impact * delta_t);
//                        
//                        double expectedNumPiecesHere        = numberOfPiecesMean[curID];      // This is the number of debris from this debIX in the catalog
//                        double probOfNoStrikeFromCurID      = pow(1. - probOfSingleStrike, expectedNumPiecesHere);
//                        double probOfNoCasualtyFromCurID    = pow(1. - probOfCasualty, expectedNumPiecesHere);
//                        double probOfNoCatastropheFromCurID = pow(1. - probOfCatastrophe, expectedNumPiecesHere);
//                        
//                        probNoStrike        *= probOfNoStrikeFromCurID;
//                        probNoCasualty      *= probOfNoCasualtyFromCurID;
//                        probNoCatastrophe   *= probOfNoCatastropheFromCurID;
//                        
//                    }
//                    
////                    double probOfAirplaneInCell = ACDensityMap[0][xindex][yindex]/numZBins; // Currently only using 2d density, so assume uniform prob over altitude
//                    
////                    double probStrike = probOfAirplaneInCell*(1. - probNoStrike)*pFail;
//                    double probStrike = (1. - probNoStrike)*pFail;
//                    // =========== End of the ACTA method ===========================================
//                    
//                    // Save the probability for this grid cell.  Only put it in the leading debIX.
//                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probImpact         = probStrike;
//                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCasualty       = (1. - probNoCasualty)*pFail;
//                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCatastrophe    = (1. - probNoCatastrophe)*pFail;
//                    
//                } } } }
//    
//    // Little bit of error checking to make sure things never went terribly awry.
//    if (minProbSingleStrike < 0.){
//        cout << "ERROR: minProbSingleStrike = " << minProbSingleStrike << endl;
//        exit(-90);  }
//    if (maxProbSingleStrike > 1.){
//        cout << "ERROR: maxProbSingleStrike = " << maxProbSingleStrike << endl;
//        exit(-90);  }
//    if (maxSweptVolume > cellVolume){
//        cout << "ERROR: maxSweptVolume = " << maxSweptVolume << endl;
//        exit(-90);  }
//    
//    
//    return;
//}


























//// ============ CUMULATIVE!!!!!!! =============
//// At some point need to take into account the curvature of the Earth...the sorting that we do in here
////   relies on the grid cell areas all being equal...so we need a map projection scheme that is equal area.
//double SkyGrid::generateAllPointsFromProbability(vector<int> numberOfPiecesMean, vector<double> arefMeanList, double thresh, double pFail){
//    //arefMeanList is no longer used.  Keep track of area stats in each cell
//    //numberOfPiecesMean is still important because it tells us how many pieces to are in each coefficient group
//    
//    double ExpectedValueCollisionSum = 0.;  // EV of collision integrated over all timesteps for the given tfail
//    
//    // Nhere
//    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh
//
//    // Clear this before we load it up
//    all_points_total.clear();
//    all_points_total.assign(NumTimestepsHere, vector<Point>());
//    
//    // Iterators for the probability grid
//    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
//
//    int tfail = 0;
//    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
//        
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//        if (it_time==ProbabilityMapDebIX.begin()){ tfail = tx; }
//        
//        vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
//        ProbabilityHere = generateProbabilityOfImpact(tx, numberOfPiecesMean, thresh, pFail);
//        
//        // Since the grid cell areas are the same AND assuming uniform aircraft density, sorting the CDF gives same ranking as sorting the EsubC values
//        // using myfunction: compares the third index (the probability) of each vector passed to it
//        // these should be sorted from highest to lowest probability
//        std::sort(ProbabilityHere.begin(), ProbabilityHere.end(), compareXYZ);
//        
//        // ~~~~~~~~ Expected Casualty Calculation ~~~~~~~~~~~~~
//        int stopIX = ProbabilityHere.size();      // index in reverse (lowest prob first)
//        double probSum = 0.;            // EsubC = probSum * ProbOfAirplane * ExpectedValueNumDebrisHere
//        //            double probLevel = thresh / (Ntotal * NumTimestepsHere * probOfAirplane * ExpectedValueNumDebrisHere * pFail);
//        double probLevel = thresh / (NumTimestepsHere - tfail + 1);     // SHOULD SUBTRACT DENOMINATOR BY TFAIL (first tx)!!!!
//
//        // Figure out where the last index is
//        while (probSum < probLevel){
//            stopIX--;                   // Subtract 1 before accessing vector for first time
//            probSum += ProbabilityHere[stopIX][3];
//            if (stopIX == 0){ break; }
//            //                cout << "probSum[" << stopIX << "] = " << probSum << "  ";
//        }
//
//        // Two options here that make stopIX == 0:
//        //  1.) because there are actually no points that need to be excluded
//        //  2.) the very last point (perhaps there is only one point total) puts you over the threshold
//    
//        // If that last index put us over the limit, add one back.
//        if (probSum > probLevel){
//            probSum -= ProbabilityHere[stopIX][3];
//            stopIX++;
//        }
//        
//        double ExpectedValueCollision = probSum;
//        ExpectedValueCollisionSum += ExpectedValueCollision;
////        cout << tx <<  ",  probLevel - " << probLevel << " stopIX  = " << stopIX << "  probSum = " << probSum << "  cumSum = " << ExpectedValueCollisionSum << endl;
////        cout << "Expected value of collision at this tstep = " << ExpectedValueCollision << endl;
////        cout << "\n\n\n\n";
//    
//        
//        // DEBUGGING
//        ProbabilityTotalStorage[tx] = ProbabilityHere;
//        stopIXStorage[tx] = stopIX;
//        
//        // resize the vector and return it
//        ProbabilityHere.resize(stopIX);
//        
//        // Save the points in the all_points vector
//        loadRemainingPointsIntoAllPoints(tx, ProbabilityHere);
//    }
//    
//    return ExpectedValueCollisionSum;
//}








/*! This replicates, as closely as reasonably possible, the way the FAA creates SUAs from probability of 
 *      impact debris clouds.  The
 *  TODO: Finish removing the coarsen grid stuff. All references to newDeltaXY_in, newDeltaZ_in
 */
double SkyGrid::generateAllPoints_CumulativeFAA(double thresh, int whichProb, double pFail){
    
    // Make sure that the probabilities have been generated, otherwise we cannot proceed
    if (not hazardProbabilitiesGenerated){
        cout << "ERROR!!!!  You cannot call generateAllPoints_CumulativeFAA before generateHazardProbabilities!" << endl;
        exit(-100); }
        
    // Since this value, for a given converged P_impact distribution, will vary dramatically as bin sizes are changed,
    //      it's really a pretty useless piece of information to return.  Will divide by the correction coefficient
    //      to scale it to be the "maxAllowableProbability per 4nm^2" value
    double maxAllowableProbability = 0.;

    double xref = XREF;
    double yref = YREF;
    double zref = ZREF;
    
    // Nhere
    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh
    
    // Clear this before we load it up
    all_points_total.clear();
    all_points_total.assign(NumTimestepsHere, vector<Point>());
    
    // This generates the spatial probability over the grid originally specified
    generateSpatialProbability(whichProb);
    
//    // The FAA uses a rather coarse grid, so convert to their coarse grid
//    //      These are the parameters of the coarsened grid
//    double newDeltaXY   = newDeltaXY_in;      //[km]
//    double newDeltaZ    = newDeltaZ_in;       //[km]  20km...This is higher than NASkm, but I need the values to nest, so hopefully this is fine

    //    projectSpatialProbabilityFAA(newDeltaXY, newDeltaZ);

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
    
//    int numZBins = (INTxx) ceil((NASkm - 0)/newDeltaZ);         // This should be 1
//    double probAircraftPresentInCell = uniformProbabilityValue;
//    double assumedArea = 13.72;     //4nm^2 in km^2
//    double coeff = newDeltaXY*newDeltaXY / (assumedArea * numZBins);
    
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
//                double curProb = SpatialProbabilty_Coarse[zindex][xindex][yindex] / probAircraftPresentInCell;
                double curProb = SpatialProbabilty[zindex][xindex][yindex] * pFail;
                
                if (curProb > (thresh)) {
                    
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
    
    
    return (maxAllowableProbability);
//    return (maxAllowableProbability/coeff);
}



















void SkyGrid::loadRemainingPointsIntoAllPoints(int tx, vector<vector<double> > ProbabilityHere){
    // Turn it into a points vector and return it
    double eps = 1e-4;  // Add a little bit to make sure this gets binned properly later
    
    unsigned long stopIX = ProbabilityHere.size();
    for (unsigned long ix = 0; ix < stopIX; ix++){
        Point hollowPoint;
        
        for (int dx = 0; dx < 2; dx++){
            for (int dy = 0; dy < 2; dy++){
                // BIG KLUDGE!!!  Hopefully fixes the fact that the swinging arm doesn't sort left-to-right
                double randx = (rand() % 1000)/1000000.;
                double randy = (rand() % 1000)/1000000.;
                
                double curx = ProbabilityHere[ix][0] + dx*xBinLength + randx;
                double cury = ProbabilityHere[ix][1] + dy*yBinLength + randy;
                double curz = ProbabilityHere[ix][2] + eps;
                
                hollowPoint.set_xyz(curx, cury, curz);
                all_points_total[tx].push_back(hollowPoint);
                
            } } }
}

/*
 * Currently this function calculates and returns P_{I}(i \mid t \leq J, f) where I corresponds to whichProb
 * TODO: Generalize this function to P_{I}(i \mid t \leq j, f), i.e. calculate only up to time j
 *          Probably a good idea to make f an input when you do that, so you can get f=t case easily
 *
 */
void SkyGrid::generateSpatialProbability(int whichProb/*, int j, int f*/){
    // whichProb selects between impact, casualty, and catastrophe
    // j is the INDEX of the maximum time to consider
    // f is the INDEX of the fail time that this probability corresponds to
    
    // Iterators for the probability grid
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    double maxVal = 0;
    double thisProb = 0;
    // Have to run through all of the timesteps first to figure out the total probabilities at each timestep
    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
        
        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
            int zindex = it_z->first;
            
            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                int xindex = it_x->first;
                
                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                    int yindex = it_y->first;
                    
                    // Options are PROB_IMPACT, PROB_CASUALTY, PROB_CATASTROPHE
                    if (whichProb == PROB_IMPACT){
                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][STORE_IX].probNoImpact;}
                    else if (whichProb == PROB_CASUALTY){
                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][STORE_IX].probNoCasualty;}
                    else if (whichProb == PROB_CATASTROPHE){
                        thisProb = ProbabilityMapDebIX[tx][zindex][xindex][yindex][STORE_IX].probNoCatastrophe;}
                    else {
                        cout << "You chose a bad probability option.  Exiting\n";
                        exit(-99);
                    }
                    
                    // This was probably incorrect in previous version, was += instead of *=
                    // Now though, thisProb is probOfNoEvent over all debris but only at a single point,
                    //  here we multiply all the time points together to get the probability of no event
                    //  over all times.
//                    SpatialProbabilty[zindex][xindex][yindex] *= thisProb;
                    
                    // Checking to make sure that default values of the map are zero
//                    if (SpatialProbabilty[zindex][xindex].count(yindex) == 0){
//                        printf("[%d][%d][%d][%d] = %E\n", tx, zindex, xindex, yindex, SpatialProbabilty[zindex][xindex][yindex]);
//                    }
                    
                    // TODO: Check to see if there are roundoff errors or something doing it this way versus above.
                    double prevVal = SpatialProbabilty[zindex][xindex][yindex];
                    SpatialProbabilty[zindex][xindex][yindex] = 1. - (prevVal+1.) * thisProb;
                    
                    maxVal = std::max(maxVal, SpatialProbabilty[zindex][xindex][yindex]);
                    
                } } }   // Ends loop back to it_z
    }
    
//    // Check on the probabilities
//    { // Changing scope so I can redefine the iterators
//        map<int, map<int, map<int,double> > >::iterator it_z;
//        map<int, map<int,double> >::iterator it_x;
//        map<int, double>::iterator it_y;
//        
//        for (it_z = SpatialProbabilty.begin(); it_z != SpatialProbabilty.end(); ++it_z){
//            int zindex = it_z->first;
//            for (it_x = SpatialProbabilty[zindex].begin(); it_x != SpatialProbabilty[zindex].end(); ++it_x){
//                int xindex = it_x->first;
//                for (it_y = SpatialProbabilty[zindex][xindex].begin(); it_y != SpatialProbabilty[zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//                    
//                    if (zindex < 0){
//                        printf("SpatialProbabilty[%d][%d][%d] = %E\n",
//                               zindex, xindex, yindex, SpatialProbabilty[zindex][xindex][yindex]);
//                    }
//                }}}
//    }

//    cout << "C++ maxVal = " << maxVal << endl;
    
//    cout << "\n\n\nDEBUG!!!" << endl;
//    projectSpatialProbabilityFAA();

    return;
}

 

map<int, map<int, map<int,double> > >SkyGrid::getSpatialProbabilty(){
    return SpatialProbabilty;
}


/*! Take SpatialProbability, which is defined on a relatively fine grid, and transfer it
 *      to a much coarser grid (3.5km on a side) which approximates the grid spacing the
 *      FAA uses (2n.m. on a side ~ 3.7km) when doing their risk calculations.
 *
 *      NOTE: This function updates the xBinLenght, yBinLength, and zBinHeight values!!!
 */
map<int, map<int, map<int,double> > > SkyGrid::projectSpatialProbabilityFAA(double newDeltaXY, double newDeltaZ){
    
//    // These are the parameters of the coarsened grid
//    double newDeltaXY   = 3.5;      //[km]
//    double newDeltaZ    = 20.;      //[km]  This is higher than NASkm, but I need the values to nest, so hopefully this is fine
    
    // This is where we'll store the new spatial probability for now
    map<int, map<int, map<int,double> > > SpatialProbabilty_Coarse;
    
    // Do some checking to make sure they nest with existing
    bool xNests = (fmod(newDeltaXY,xBinLength) < 1e-12);
    bool yNests = (fmod(newDeltaXY,yBinLength) < 1e-12);
    bool zNests = (fmod(newDeltaZ ,zBinHeight) < 1e-12);
    
    if (not (xNests and yNests and zNests)){
        cout << "ERROR!!!  projectSpatialProbabilityFAA does not nest\n";
        printf("%f, %f, %f\n", fmod(newDeltaXY,xBinLength), fmod(newDeltaXY,yBinLength), fmod(newDeltaZ ,zBinHeight));
        return SpatialProbabilty_Coarse;
    }

    // These are the iterators for map<int, map<int, map<int,double> > > SpatialProbabilty
    map<int, map<int, map<int,double> > >::iterator it_z;
    map<int, map<int,double> >::iterator it_x;
    map<int,double>::iterator it_y;
    
    for (it_z = SpatialProbabilty.begin(); it_z != SpatialProbabilty.end(); ++it_z){
        int zindex      = it_z->first;   // This index may need to be coarsened since the FAA doesn't use different altitudes
        int newZIndex   = floor((ZREF + zindex * zBinHeight)/newDeltaZ);
        
        // printf("zindex = %d, newZIndex = %d, zNests = %s\n", zindex, newZIndex, zNests ? "True":"False");
        
        for (it_x = it_z->second.begin(); it_x != it_z->second.end(); ++it_x){
            int xindex      = it_x->first;
            int newXIndex   = floor((XREF + xindex*xBinLength)/newDeltaXY);
            
            // printf("    xindex = %d, newXIndex = %d, xNests = %s\n", xindex, newXIndex, xNests ? "True":"False");

            
            for (it_y = it_x->second.begin(); it_y != it_x->second.end(); ++it_y){
                int yindex      = it_y->first;
                int newYIndex   = floor((YREF + yindex*yBinLength)/newDeltaXY);
                
                SpatialProbabilty_Coarse[newZIndex][newXIndex][newYIndex] += it_y->second;
                // printf("        yindex = %d, newYIndex = %d, yNests = %s\n", yindex, newYIndex, yNests ? "True":"False");

            }
            
        }
        
        
    }
    
    // Update the values
    SpatialProbabilty = SpatialProbabilty_Coarse;
    xBinLength = newDeltaXY;
    yBinLength = newDeltaXY;
    zBinHeight = newDeltaZ;
    
    return SpatialProbabilty_Coarse;
    
    
}








void SkyGrid::DumpGridToMatlab(char *fileName){
    
    // Check that you can still use it
    if (!isProbability){
        cout << "ERROR, youre trying to use DumpGridToMatlab before converting to probabilities\n\n";
        exit(-13);
    }
    
    int debugTime = 40 ;
//    int debugTime = 60 ;
    
    ofstream outfile;
	outfile.open(fileName, ios::out);
    
    // Set the precision and make sure output is set to scientific
	outfile.precision(9);
	outfile << std::scientific;
    
    // Iterators for the probability grid
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    outfile << "clear all; close all; clc;\n\n";
    
//    cout << "DUMPING ProbabilityMapDebIX\n";
    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
        
        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
        
//        if (tx != debugTime) {continue; }
        
        outfile << "vec{" << tx << "} = [\n";

        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
            int zindex = it_z->first;
            
            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                int xindex = it_x->first;
                
                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                    int yindex = it_y->first;
                    
                    double probSumHere = 0.;
                    for (it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin(); it_ID != ProbabilityMapDebIX[tx][zindex][xindex][yindex].end(); ++it_ID){
                        int curID = it_ID->first;
                        
                        binData PD = it_ID->second;
                        
                        outfile << zindex << "  " << xindex << "  " << yindex << "  " << curID <<  "  " << PD.probDebris << "  " << PD.avgVel
                            << "  " << PD.minVel << "  " << PD.maxVel << ";\n";
                        
//                        double probOfStrike = (arefMeanList[curID]*m2tokm2 + airplaneArea)/cellArea;
//                        probSumHere += (it_ID->second) * numberOfPiecesMean[curID] * probOfStrike;
                        
                        
                    } } } }
        outfile << "];\n\n\n\n";
    }
    
    
    
//    temp4Vec[0] = lowerleftX;
//    temp4Vec[1] = lowerleftY;
//    temp4Vec[2] = lowerleftZ;
//    temp4Vec[3] = 1. - probNoStrike;
    
    map<int, vector<vector<double> > >::iterator it_StorageTime;
//    cout << "DUMPING ProbabilityTotalStorage\n";

//    int tsteps = ProbabilityTotalStorage.size();
//    for (int tx = 0; tx < tsteps; tx++){
    for (it_StorageTime = ProbabilityTotalStorage.begin(); it_StorageTime != ProbabilityTotalStorage.end(); ++it_StorageTime){
        int tx = it_StorageTime->first;
//        cout << "tx" << endl;
        
//        if (tx != debugTime) {continue; }
        
        outfile << "ProbabilityTotalStorage{" << tx << "} = [\n";

        int numHere = ProbabilityTotalStorage[tx].size();
        for (int ix = 0; ix < numHere; ix++){
            
            // Outputing in the format z,x,y,prob
            outfile << ProbabilityTotalStorage[tx][ix][2] << "  " << ProbabilityTotalStorage[tx][ix][0] << "  "
                << ProbabilityTotalStorage[tx][ix][1] << "  " << ProbabilityTotalStorage[tx][ix][3] << ";\n";

        }
        outfile << "];\n\n\n\n";
    }
    
//    cout << "DUMPING stopIXStorage\n";
    map<int, int>::iterator it_StopIXTime;

    for (it_StopIXTime = stopIXStorage.begin(); it_StopIXTime != stopIXStorage.end(); ++it_StopIXTime){
        int tx = it_StopIXTime->first;
        
//        if (tx != debugTime) {continue; }
        
        outfile << "stopIXStorage{" << tx << "} = " << stopIXStorage[tx] << ";\n";

    }
    outfile << "\n\n";
    
    outfile <<  "NumTimestepsHere = " << getNumRange() << endl;

    outfile.close();
    
    return;
}





void SkyGrid::GoMatlab(string fileName, vector<Point> tempVec){
    
    // Check that you can still use it
    if (isProbability){
        cout << "ERROR, youre trying to use GoMatlab after converting to probabilities\n\n";
        exit(-13);
    }
    
    ofstream outfile;
	outfile.open(fileName.c_str(), ios::out);
    
    // Set the precision and make sure output is set to scientific
	outfile.precision(9);
	outfile << std::scientific;
    
    outfile << "clear all; close all; clc;\n\n";
    
    outfile << "vec = [\n";
    //    for (int tx = (all_points_num_range-1); tx < all_points_num_range; tx ++){
    //        int num_current_points = (INTxx) pts[tx][curIX].size();
    for (int ix = 0; ix < tempVec.size(); ix++){
        outfile << tempVec[ix].get_x() << "  " << tempVec[ix].get_y() << ";\n";
    }
    //}
    
    outfile << "];\n\n\n\n";
    outfile << "scatter(vec(:,1), vec(:,2))\n";
    
    outfile.close();    
    
    return;
}




int SkyGrid::identifyYourself(){
//    cout << "Derived class SkyGrid returns 1" << endl;
    return 1;
}

double SkyGrid::getZBinHeight(){
    return zBinHeight;
}

int SkyGrid::getNumRange(){
    int ans = -1;
    if (isProbability){
//        ans = Probability.size();
        
        // Check if the vector is empty
        if (ProbabilityMapDebIX.begin() == ProbabilityMapDebIX.end()){
            ans = 0;
        } else {
            int firstIX = std::min(0, ProbabilityMapDebIX.begin()->first);  // the index of the least timestep
            int lastIX = (--(ProbabilityMapDebIX.end()))->first;  // the index of the greatest timestep
            ans = lastIX + 1 + abs(firstIX);    // firstIX will be <= 0

        }
    
//        // Note that end() puts us PAST the end of the container.  We must decrement once to be at the final element
//        cout << "Segfault right here\n";
//        int firstIX = std::min(0, ProbabilityMapDebIX.begin()->first);  // the index of the least timestep
//
//        cout << "truth = " << (ProbabilityMapDebIX.begin() == ProbabilityMapDebIX.end() )  << endl;
//        
////        cout << "ProbabilityMapDebIX.begin() = " << *(ProbabilityMapDebIX.begin())  << endl;
////        cout << "ProbabilityMapDebIX.end() = " <<   *(ProbabilityMapDebIX.end())    << endl;
////        cout << "--ProbabilityMapDebIX.end() = " << --(ProbabilityMapDebIX.end())   << endl;
//        cout << "Segfault right here\n";
//        
//        
//        int lastIX = (--(ProbabilityMapDebIX.end()))->first;  // the index of the greatest timestep
//
//        ans = lastIX + 1 + abs(firstIX);    // firstIX will be <= 0
        
        
//        cout << "ISPRObABLITY ========================================================= " << ans << "   " << lastIX;
//        cout << "\n\n\n\n\n\n\n\n";
    }
    else {
//        ans = Grid.size();
//        cout << "In moving to DebIX versions of things, this here might be broken.  Exiting." << endl;
//        exit(17);
        
        int lastIX = (--totalNumPtsAtStepMapDebIX.end())->first;  // the index of the greatest timestep
        int firstIX = std::min(0, totalNumPtsAtStepMapDebIX.begin()->first);  // the index of the least timestep
        ans = lastIX - firstIX + 1; // this is fine because firstIX will be <= 0
    }
    
    return ans;
}


// Eqn 5.7
double biweightKernel(double t){
    return ( 15./16. * fmax( std::pow(1 - t*t,2), 0.)  );
}




/*!
 * Doxygen looks like this
 */
vector<map<int,SkyGrid::binData> > SkyGrid::ASHDesiredPoint(double h1_in, double h2_in, vector<vector<double> > desiredPts){

    
    // ==== Some ASH setup stuff ====
    
    // h is the smoothing parametere applied in each direction
    double h1 = h1_in;    // Assuming that 10km is a natural spreading distance
    double h2 = h2_in;
    
    double delta1 = xBinLength;     // Size of initial histogram cells
    double delta2 = yBinLength;
    
    int m1 = ceil(h1/delta1);     // Number of histograms to average over in X direction
    int m2 = ceil(h2/delta2);     // Number of histograms to average over in Y direction
    
    // Update the spreading distance to be consistent
    h1 = m1 * delta1;
    h2 = m2 * delta2;
    
    // Allocate the weights assuming a bivariate kernel
    map<int, double> w_m1;          // The actual weights
    map<int, double> w_m2;
    
    // Find the normalization constants and then initialize the weight vectors
    double normalizationConst1 = 0.;
    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
        normalizationConst1 += biweightKernel(ix/m1); }
    
    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
        w_m1[ix] = m1 * delta1 * biweightKernel(ix/m1)/normalizationConst1; }
    
    double normalizationConst2 = 0.;
    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
        normalizationConst2 += biweightKernel(ix/m2); }
    
    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
        w_m2[ix] = m2 * delta2 * biweightKernel(ix/m2)/normalizationConst2; }
    
//    double normFactor = ((double) numDebrisPerIXSimulated) * h1 * h2;
    
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_ID;
    map<int, map<int, map<int,binData> > >::iterator it_z;
    map<int, map<int, binData> >::iterator it_x;
    map<int, binData>::iterator it_y;
    
    
    // ==== Now jump into finding the probabilities ====
    
    // Loop over the points
    int numPts = desiredPts.size();
    
    // Return a map of probability datas (for each single point)
    vector<map<int,binData> > probabilityRecord;
    probabilityRecord.assign(numPts, map<int,binData>());
    
    for (int ix = 0; ix < numPts; ix++){
    
        // Unpack the desired point
        int tx0         = floor((desiredPts[ix][0] - 0)/getDeltaT());
        int zindex0     = floor((desiredPts[ix][1] - ZREF)/zBinHeight);
        int xindex0     = floor((desiredPts[ix][2] - XREF)/xBinLength);
        int yindex0     = floor((desiredPts[ix][3] - YREF)/yBinLength);
        
        // Look for iterator to desired timestep
        it_time=GridMapDebIX.find(tx0);
        
        // Does the time exist?
        if (it_time != GridMapDebIX.end()){
            // Yes, it's in the grid.
            
            // I recall having a good reason for organizing the grid map this way, but right now it seems pretty stupid to
            //      have the debrisID be right here.  Must loop through them all.
            
            for (it_ID = it_time->second.begin(); it_ID != it_time->second.end(); ++it_ID){
                
                //For this ID, look for a z-iterator
                it_z = it_ID->second.find(zindex0);
                
                // Does it exist?
                if (it_z != it_ID->second.end()){
                    
                    
    //                double avgVel = 0.;
    //                double minVel = 0.;
    //                double maxVel = 0.;
    //
    //                double avgMass = 0.;
    //                double minMass = 0.;
    //                double maxMass = 0.;
    //                
    //                double avgArea = 0.;
    //                double minArea = 0.;
    //                double maxArea = 0.;
    //                
    //                int numTouched = 0;
                    
                    binData PD;
                    
                    int curID = it_ID->first;
                    double normFactor = ((double) totalNumPointsPassedInPerID[curID]) * h1 * h2;
                    
                    // Yes, it does.  Now need to look for iterators to all of the x,y locations that the ASH would feed into this point
                    for (int i = 1 - m1; i <= m1 - 1; i++){
                        // Note that i starts negative (left) and scans positive (right)
                        
                        // Look for the x-iterator
                        it_x = it_z->second.find(xindex0 + i);
                        
                        // Does it exist?
                        if (it_x != it_z->second.end()){
                            
                            
                            // Yes it does.  Now loop through the possible y iterators
                            for (int j = 1 - m2; j <= m2 - 1; j++){
                                
                                it_y = it_x->second.find(yindex0 + j);
                                
                                // Does it exist?
                                if (it_y != it_x->second.end()){
                                    
                                    //YES!  So calculate this point's contribution to the probability at the desired point
                                    // By virtue of this iterator existing, that means there is a non-zero histogram here
                                    binData gridData = it_y->second;
                                    
                                    // Spread the velocity around
                                    PD.avgVel = (PD.numTouched*PD.avgVel + gridData.probDebris*gridData.avgVel)/(PD.numTouched + gridData.probDebris);
                                    PD.maxVel = std::max(PD.maxVel, gridData.maxVel);
                                    PD.minVel = std::min(PD.minVel, gridData.minVel);
                                    
                                    // Spread the mass around
                                    PD.avgMass = (PD.numTouched*PD.avgMass + gridData.probDebris*gridData.avgMass)/(PD.numTouched + gridData.probDebris);
                                    PD.maxMass = std::max(PD.maxMass, gridData.maxMass);
                                    PD.minMass = std::min(PD.minMass, gridData.minMass);
                                    
                                    // Spread the area around
                                    PD.avgArea = (PD.numTouched*PD.avgArea + gridData.probDebris*gridData.avgArea)/(PD.numTouched + gridData.probDebris);
                                    PD.maxArea = std::max(PD.maxArea, gridData.maxArea);
                                    PD.minArea = std::min(PD.minArea, gridData.minArea);
                                    
                                    PD.numTouched += gridData.probDebris;
                                    
                                    // Do the actual ASH with the probability (indexes of the weights will be negative...mirror image kinda)
                                    PD.probDebris += gridData.probDebris * w_m1[-i] * w_m2[-j] / normFactor;
                                    
                                } } } }
                        
                    // Create this value in the map and store it
                    probabilityRecord[ix][it_ID->first] = PD;
                } } }
    }
    
    
    
    
    return probabilityRecord;
}


void SkyGrid::ASH2(double h1_in, double h2_in){
    // Bivariate Average Shifted Histogram
    // Based on the book Multivariate Density Estimation by DAVID W. SCOTT
    // Scott chooses the weights to get a discrete approximation of the density function,
    //      but I'm choosing them so that I get a mass function.  Could just delete a lot of the normalization stuff actually
    
    // I think the way this should ACTUALLY work is that you pick the h's and then the m's
    //      derive from that.  This algo does it the other way around.  CHANGE THIS EVENTUALLY.
    
    // Check that you can still use it
    if (isProbability){
        cout << "ERROR, youre trying to use ASH2 after converting to probabilities\nDo Nothing\n";
        exit(-16);
        return;
    }
    
    isProbability = true;
    
    // Since everything here is anticpated to either work or blow up...
    doneASH = true;
    
    // h is the smoothing parametere applied in each direction
    double h1 = h1_in;    // Assuming that 10km is a natural spreading distance
    double h2 = h2_in;
    
    double delta1 = xBinLength;     // Size of initial histogram cells
    double delta2 = yBinLength;
    
    int m1 = ceil(h1/delta1);     // Number of histograms to average over in X direction
    int m2 = ceil(h2/delta2);     // Number of histograms to average over in Y direction
    
    // Update the spreading distance to be consistent
    h1 = m1 * delta1;
    h2 = m2 * delta2;
    
    //    double delta1 = xBinLength;     // Size of initial histogram cells
    //    double delta2 = yBinLength;
    //
    //    int m1 = 5;     // Number of histograms to average over in X direction
    //    int m2 = 5;     // Number of histograms to average over in Y direction
    //
    //    double h1 = m1 * delta1;
    //    double h2 = m2 * delta2;
    
    // Allocate the weights assuming a bivariate kernel
    map<int, double> w_m1;          // The actual weights
    map<int, double> w_m2;
    
    // Find the normalization constants and then initialize the weight vectors
    double normalizationConst1 = 0.;
    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
        normalizationConst1 += biweightKernel(ix/m1); }
    
    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
        w_m1[ix] = m1 * delta1 * biweightKernel(ix/m1)/normalizationConst1; }
    
    double normalizationConst2 = 0.;
    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
        normalizationConst2 += biweightKernel(ix/m2); }
    
    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
        w_m2[ix] = m2 * delta2 * biweightKernel(ix/m2)/normalizationConst2; }
    
    // This is to explicitly bound the probabilities by zeros and is my own addition in order
    //      to make the map structure work the way i want for interpolating values later on.
    int extraTouch = 1;     // not set up to be anything other than 0 or 1
    // Actually, i get different risk numbers if i use extraTouch at all.  Why?  Don't use until fixed.
    // I assume it's because of the spreading around velocities, areas, etc.
    // You'll create zero-probability areas that have non-zero velocities and whatnot...is that actually a problem???
    // YES!  Imagine your extraTouching a cell that already has had a probability assigned from somewhere else, then you'll
    //      alter all of the other statistics EXCEPT the probability...which will change the risk answers.
    if (extraTouch == 1){
        w_m1[1 - m1 - extraTouch] = 0.;
        w_m1[m1 - 1 + extraTouch] = 0.;
        w_m2[1 - m2 - extraTouch] = 0.;
        w_m2[m2 - 1 + extraTouch] = 0.;
    }
    
    //    map<int, map<int, map<int, map<int,double> > > > ProbabilityMap;
    //    map<int, map<int,double> > SubProbabilityMap;
    
    {   // Starting a new scope so that the iterators will fall out of scope and i can use the same names again in a minute
        map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
        map<int, map<int, map<int, map<int,binData> > > >::iterator it_ID;
        map<int, map<int, map<int,binData> > >::iterator it_z;
        map<int, map<int, binData> >::iterator it_x;
        map<int, binData>::iterator it_y;
        
        //    cout << "Start ASH...";
        // Spread the probabilities around with the weights
        for (it_time=GridMapDebIX.begin(); it_time!=GridMapDebIX.end(); ++it_time) {
            int tx = it_time->first;
            
//            cout << "ash2 tx = " << tx << "/" << GridMapDebIX.size() << endl;
            
            for (it_ID = GridMapDebIX[tx].begin(); it_ID != GridMapDebIX[tx].end(); ++it_ID){
                int curID = it_ID->first;
                
                double normFactor = ((double) totalNumPointsPassedInPerID[curID]) * h1 * h2;
                
                for (it_z = GridMapDebIX[tx][curID].begin(); it_z != GridMapDebIX[tx][curID].end(); ++it_z){
                    int zindex = it_z->first;
                    
                    for (it_x = GridMapDebIX[tx][curID][zindex].begin(); it_x != GridMapDebIX[tx][curID][zindex].end(); ++it_x){
                        int xindex = it_x->first;
                        
                        for (it_y = GridMapDebIX[tx][curID][zindex][xindex].begin(); it_y != GridMapDebIX[tx][curID][zindex][xindex].end(); ++it_y){
                            int yindex = it_y->first;
                            // By virtue of this iterator existing, that means there is a non-zero histogram here
                            binData gridData = GridMapDebIX[tx][curID][zindex][xindex][yindex];
                            
                            // Loop over all the weights and spread the probabilities around
                            for (int i = 1 - m1 - extraTouch; i <= m1 - 1 + extraTouch; i++){
                                for (int j = 1 - m2 - extraTouch; j <= m2 - 1 + extraTouch; j++){
                                    // NOTE: The index order of the probability map is different from that of the grid map!!!!
                                    binData PD = ProbabilityMapDebIX[tx][zindex][xindex + i][yindex + j][curID];
                                    
                                    // Do nothing if we're extraTouching.  Just creating the map cell in the line above is enough.
                                    if ((w_m1[i] == 0.) or (w_m2[j] == 0.)){
                                        // Do Nothing
                                    } else {
                                        // Spread the velocity around
                                        PD.avgVel = (PD.numTouched*PD.avgVel + gridData.probDebris*gridData.avgVel)/(PD.numTouched + gridData.probDebris);
                                        PD.maxVel = std::max(PD.maxVel, gridData.maxVel);
                                        PD.minVel = std::min(PD.minVel, gridData.minVel);
                                        
                                        // Spread the mass around
                                        PD.avgMass = (PD.numTouched*PD.avgMass + gridData.probDebris*gridData.avgMass)/(PD.numTouched + gridData.probDebris);
                                        PD.maxMass = std::max(PD.maxMass, gridData.maxMass);
                                        PD.minMass = std::min(PD.minMass, gridData.minMass);
                                        
                                        // Spread the area around
                                        PD.avgArea = (PD.numTouched*PD.avgArea + gridData.probDebris*gridData.avgArea)/(PD.numTouched + gridData.probDebris);
                                        PD.maxArea = std::max(PD.maxArea, gridData.maxArea);
                                        PD.minArea = std::min(PD.minArea, gridData.minArea);
                                        
                                        PD.numTouched += gridData.probDebris;
                                        
                                        // Do the actual ASH with the probability
                                        PD.probDebris += gridData.probDebris * w_m1[i] * w_m2[j] / normFactor;
                                        
                                        ProbabilityMapDebIX[tx][zindex][xindex + i][yindex + j][curID] = PD;

                                    }
                                    
                                    
                                    
                                }}
                        }}
                }}}
    }

    // Check the probabilities (Be safe and keep this in.  Don't be an asshole and delete such an important check)
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    for (it_time=ProbabilityMapDebIX.begin(); it_time!=ProbabilityMapDebIX.end(); it_time++) {
        int tx = it_time->first;
        
        // At every timestep, the probabilities should sum to 1
        map<int, double> checkSum;
        
        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
            int zindex = it_z->first;
            
            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                int xindex = it_x->first;
                
                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                    int yindex = it_y->first;
                    
                    for (it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin(); it_ID != ProbabilityMapDebIX[tx][zindex][xindex][yindex].end(); ++it_ID){
                        int curID = it_ID->first;
                        
                        checkSum[curID] += ProbabilityMapDebIX[tx][zindex][xindex][yindex][curID].probDebris;
                    } } } }
        
        map<int, double>::iterator it_CS;
        for (it_CS = checkSum.begin(); it_CS != checkSum.end(); ++it_CS){
            int curID = it_CS->first;

            // printf("CheckSum[%d][%d] = %E \n", tx, curID, checkSum[curID]-1.);
            if (fabs(checkSum[curID]-1.) > 1e-13){
                printf("CheckSum[%d][%d] = %E \n", tx, curID, checkSum[curID]-1.);
                exit(-15);}

        }
        
    }
    
    return;
}









// This function is mainly only for debugging purposes
map<double, map<double, map<double,double> > >SkyGrid::SendGridToPython(int tx_desired){
    
//    cout << "Entering SendGridToPython" << endl;

    map<double, map<double, map<double,double> > > GridAsVector;   //[z][x][y]
//    GridAsVector[-1.][-2.][-3.] = 666.;
    
    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh
//    cout << "NumTimestepsHere = " << NumTimestepsHere << endl;
    
    // Clear this before we load it up
//    GridAsVector.assign(NumTimestepsHere, vector<double>());
    
//    int numDebrisCategories = countDebrisIX.size();
    int numDebrisCategories = numPtsLessThanReactionInNAS.size();
    
    // First must dump the debIX grids down to a single grid
    map<int, map<int, map<int, map<int,double> > > >CollapsedMap;
    
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    // Find the probability vector at this timestep for all points in the grid.
    //    vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
    
    //    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
    
    Point tempPt;
    
    it_time = ProbabilityMapDebIX.find(tx_desired);
    if (it_time != ProbabilityMapDebIX.end()){
        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
        
        
//        // For this area, we need to find the min and max latitudes
//        // These are the indicies which are equally spaced
//        int minx = 1e9;
//        int maxx = -1e9;
//        int miny = 1e9;
//        int maxy = -1e9;
        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
            int zindex = it_z->first;
//
//            // For each z-level, find the first (min) xval
//            it_x = ProbabilityMapDebIX[tx][zindex].begin();
//            minx = std::min(it_x->first, minx);
//            
//            // For each z-level, find the last (max) xval
//            it_x = --(ProbabilityMapDebIX[tx][zindex].end());
//            maxx = std::max(it_x->first, maxx);
//            
//            // For each z-level, find the first (min) yval
//            it_x = ProbabilityMapDebIX[tx][zindex].begin();
//            minx = std::min(it_x->first, minx);
//            
//            // For each z-level, find the last (max) yval
//            it_x = --(ProbabilityMapDebIX[tx][zindex].end());
//            maxx = std::max(it_x->first, maxx);
//        }
//        
//        cout << "minx,max = " << minx << ", " << maxx << endl;
            
            
        
            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
                int xindex = it_x->first;

                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
                    int yindex = it_y->first;

                    double curProb = 0.;
                    for (it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin(); it_ID != ProbabilityMapDebIX[tx][zindex][xindex][yindex].end(); ++it_ID){
                        int curID = it_ID->first;

                        binData PD = (it_ID->second);
                        double probDebrisHere = ((it_ID->second).probDebris);

                        curProb += probDebrisHere;

//                        CollapsedMap[tx][zindex][xindex][yindex] += probDebrisHere/numDebrisCategories;

                    }
                    
                    double lowerleftX = XREF + xindex*xBinLength;   // Note that xindex is most likely negative
                    double lowerleftY = YREF + yindex*yBinLength;   // Note that yindex may be positive or negative
                    double lowerleftZ = ZREF + zindex*zBinHeight;
                    tempPt.set_xyz(lowerleftX + 0.5*xBinLength, lowerleftY + 0.5*yBinLength, lowerleftZ);
                    
                    double centerLat = tempPt.get_gdLatDeg();
                    double centerLon = tempPt.get_lonDeg();

                    GridAsVector[lowerleftZ][centerLon][centerLat] = curProb;

                    // =========== End of the ACTA method ===========================================

//                    // Save the probability for this grid cell.  Only put it in the leading debIX.
//                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probImpact         = probStrike;
//                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCasualty       = probOfAirplaneInCell*(1. - probNoCasualty)*pFail;
//                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCatastrophe    = probOfAirplaneInCell*(1. - probNoCatastrophe)*pFail;

                } } }
    } else {
        cout << "This time does not exist" << endl;
    }
    
    // Then format the grid into a vector that can be passed back
    
    
//    cout << "Leaving SendGridToPython" << endl;
    return GridAsVector;
}




void SkyGrid::UploadAircraftTrackMap(map<int, pair<vector<vector<double> >, string> > AircraftTrackMap_in){
//    int len = inputVec.size();
//    cout << "len of incoming vec = " << len << endl;
//    for (int ix = 0; ix < len; ix++){
//        cout << "inputVec = " << inputVec[ix] << endl;
//    }
    
    AircraftTrackMap = AircraftTrackMap_in;
    
    // From python aircraftRecord[acid].append([float(curTrackTime), curLat, curLon, curLevel, curSpeed])
    
//    map<int, vector< vector<double> > >::iterator it_ACID;
//    for (it_ACID = AircraftTrackMap.begin(); it_ACID != AircraftTrackMap.end(); ++it_ACID){
//        int acid = it_ACID->first;
////        string acType = it_ACID->second["acType"];
////
////        cout << "acid = " << acid << " is type " << acType << endl;
//        
//        
//        int lenHere = (it_ACID->second).size();
//        for (int ix = 0; ix < std::min(lenHere,5); ix++){
//            cout << "   " << (it_ACID->second)[ix][1] << "  " << (it_ACID->second)[ix][2] << endl;
//        }
//    }
    
    
    return;
}


void SkyGrid::UploadAircraftPropertiesMap(map<string,map<string,double> > AircraftPropertiesMap_in){
    
    AircraftPropertiesMap = AircraftPropertiesMap_in;
    
//    map<string,map<string,double> >::iterator it_acType;
//    map<string,double>::iterator it_property;
//    
//    cout << "========= Dumping the AircraftPropertiesMap ========= " << endl;
//    for (it_acType = AircraftPropertiesMap.begin(); it_acType != AircraftPropertiesMap.end(); ++it_acType){
//        cout << it_acType->first << endl;
//        for (it_property = it_acType->second.begin(); it_property != it_acType->second.end(); ++it_property){
//            cout << "   " << it_property->first << " = " << it_property->second << endl;
//        }
//    }
    
    return;
}





map<int, double> SkyGrid::CalculateRiskToIndividualAircraft_OnTheFly(vector<int> numberOfPiecesMean, vector<double> arefMean, int secondsFromMidnightUTC,
                                                                     double h1_in, double h2_in){
    // NOTE: I'm going to require that the aircraft data come at a fixed timestep so if you ever wind up
    //   having an irregular timestep then this will break.  Don't break it.
    int planeStepSec = 1;
    
    double delta_t = getDeltaT();
    double ft_2_km = 0.0003048;
    
    // Open the file to dump the debug output
    string debugFileName = "debugOutput.txt";
    ofstream outfile;
	outfile.open(debugFileName.c_str(), ios::out);
    
    // Dump a little bit of the AircraftTrackMap to stdout
    map<int, pair< vector< vector<double> >, string> >::iterator it_ACID;
    for (it_ACID = AircraftTrackMap.begin(); it_ACID != AircraftTrackMap.end(); ++it_ACID){
        int acid        = it_ACID->first;           // the aircraft id -- self-explanatory
        string acType   = it_ACID->second.second;   // B737, etc
        
        cout << "acid = " << acid << ", model = " << acType << endl;
        int lenHere = (it_ACID->second.first).size();
        
        // Find the index where the explosion occurs
        int tx = 0;
        while (tx < lenHere){
            if ((it_ACID->second.first)[tx][0] == secondsFromMidnightUTC) {
                break; }
            
            tx++;
        }
        
        for (int ix = 0; (tx + ix) < std::min(lenHere,tx + 5); ix++){
            cout << "   " << (it_ACID->second.first)[tx + ix][0] << "   "
            << (it_ACID->second.first)[tx + ix][1] << "  "
            << (it_ACID->second.first)[tx + ix][2] << "  "
            << (it_ACID->second.first)[tx + ix][3] << "  "
            << (it_ACID->second.first)[tx + ix][4] << "  "
            << endl;
        }
    }
    
    // Dump a little bit of the AircraftPropertiesMap to stdout
    map<string,map<string,double> >::iterator it_acType;
    map<string,double>::iterator it_property;
    
    cout << "========= Dumping the AircraftPropertiesMap ========= " << endl;
    for (it_acType = AircraftPropertiesMap.begin(); it_acType != AircraftPropertiesMap.end(); ++it_acType){
        cout << it_acType->first << endl;
        for (it_property = it_acType->second.begin(); it_property != it_acType->second.end(); ++it_property){
            cout << "   " << it_property->first << " = " << it_property->second << endl;
        }
    }
    
    // Start the calculation
    
//    // The iterators for ProbabilityMapDebIX
//    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
//    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
//    map<int, map<int, map<int,binData> > >::iterator it_x;
//    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    double cellVolume = xBinLength*yBinLength*zBinHeight;
    
    map<int, double> probabilityOfImpactRecord;
    
    // Loop over the aircraft themselves
    // NEED THE TIME OFFSET!!!
    
    for (it_ACID = AircraftTrackMap.begin(); it_ACID != AircraftTrackMap.end(); ++it_ACID){
        int acid        = it_ACID->first;           // the aircraft id -- self-explanatory
        string acType   = it_ACID->second.second;   // B737, etc
        
        probabilityOfImpactRecord[acid] = 0.;
        
        cout << "acid = " << acid << ", model = " << acType << endl;
        int lenHere = (it_ACID->second.first).size();
        

        
        // Find the areas that are relevant for this aircraft
        //        acClass = 1
        //        frontArea = 0.0223973
        //        topArea = 0.154464
        int     acClass     = AircraftPropertiesMap[acType]["acClass"];
        double  frontArea   = AircraftPropertiesMap[acType]["frontArea"];   // These areas come in km^2
        double  topArea     = AircraftPropertiesMap[acType]["topArea"]; 
        

        
        // Find the index where the explosion occurs.  Time Explosion Offest = TEO
        // Find the index where the plane is in the air AND the explosion occurs / has previously occurred
        int plane_tx = 0;
        while (plane_tx < lenHere){
            if ((it_ACID->second.first)[plane_tx][0] >= secondsFromMidnightUTC) {
                break; }
            
            plane_tx++; }
        
        if (plane_tx == lenHere){
            // This aircraft is not in the air at the same time as the accident
            continue;
        }
        
        int TEO = (it_ACID->second.first)[plane_tx][0];
        cout << "TEO = " << TEO << ",  secondsFromMidnightUTC = " << secondsFromMidnightUTC << endl;

        double secondsSinceImpact   = 0.;

        Point tempPt;
        for (plane_tx = plane_tx; plane_tx < lenHere; plane_tx++){
            
            // What's the time index at this point?
            double plane_seconds       = (it_ACID->second.first)[plane_tx][0];
//            int debris_tx           = (int) floor((plane_seconds - TEO)/all_points_delta_t);
            int debris_tx           = (int) floor((plane_seconds - secondsFromMidnightUTC)/all_points_delta_t);
            double aircraftLat      = (it_ACID->second.first)[plane_tx][1];
            double aircraftLon      = (it_ACID->second.first)[plane_tx][2];
            double aircraftZ        = (it_ACID->second.first)[plane_tx][3];
            double aircraftSpeed    = (it_ACID->second.first)[plane_tx][4];

            // This should be the xyz coordinate of the map
            tempPt.set_xy_from_latlon(aircraftLat * PI/180., aircraftLon * PI/180.);

            //            int xindex = floor((tempPt.get_x() - XREF)/xBinLength);
            //            int yindex = floor((tempPt.get_y() - YREF)/yBinLength);
            //            int zindex = floor((aircraftZ - ZREF)/zBinHeight);

            
            // Look up the probability at this point.  Only do one point at a time, for the moment, to minimize changes to code flow.
            vector<vector<double> > ptsOfInterest;
            ptsOfInterest.assign(1,vector<double>());
            ptsOfInterest[0].assign(4,0.);
            
//            ptsOfInterest[0][0] = (plane_seconds - TEO);
            ptsOfInterest[0][0] = (plane_seconds - secondsFromMidnightUTC);
            ptsOfInterest[0][1] = aircraftZ;
            ptsOfInterest[0][2] = tempPt.get_x();
            ptsOfInterest[0][3] = tempPt.get_y();
            
            // Call the function
            vector<map<int,binData> > desiredProbabilities =
                ASHDesiredPoint(h1_in, h2_in, ptsOfInterest);

            
            // Only had one point, so no need to loop, but anticipating a loop sometime
            int numPtsOfInterest = ptsOfInterest.size();
            int curPt = 0;
            map<int,binData> curProbMap = desiredProbabilities[curPt];
            
            double probNoStrike         = 1.;
            

            // If there was no overlap between debris and aircraft at this timestep, then curProbMap should be empty and this loop be skipped
            for (it_ID = curProbMap.begin(); it_ID != curProbMap.end(); ++it_ID){
                int curID = it_ID->first;
                binData PD = (it_ID->second);
                // Assuming that if we know there is debris in this cell, the location of that debris is uniformly likely to be anywhere in the volume.
                
                double probDebrisHere = (PD.probDebris)/cellVolume;
                
                
//                // ==== Old ACTA way ====
//                double AF = pow(d_Airplane_front + sqrt(PD.avgArea),2);
//                double AT = pow(d_Airplane_top + sqrt(PD.avgArea),2);
//                double probOfSingleStrike = probDebrisHere * (AF*speed787 + AT*PD.avgVel) * delta_t;
                
                // ==== New ACTA way ====
                // Equation takes mass in grams, outputs ft^2, so need to convert to km^2
                double theta = atan2(aircraftSpeed, PD.avgVel);  // Sure hope velocity is in km/s
                double A_Proj = frontArea*sin(theta) + topArea*cos(theta);
                double V_impact = sqrt(pow(aircraftSpeed, 2) + pow(PD.avgVel,2));
                
                double A_Casualty   = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);   // This is the case for pieces over 300g
                double A_Catastrope = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);
                
                if (PD.avgMass < 0.001){
                    // Throw out pieces that are less than 1g.  They pose no danger
                    A_Casualty      = 0.;
                    A_Catastrope    = 0.;
                }
                else if (PD.avgMass < 0.300){
                    A_Casualty      = ((0.0085*pow(PD.avgMass * 1e3,2) + 8.5*(PD.avgMass * 1e3) + 200)) * pow(ft_2_km,2);
                    A_Catastrope    = (0.025 * pow(PD.avgMass * 1e3,2) + 4*(PD.avgMass * 1e3)) * pow(ft_2_km,2);
                }
                
                double probOfSingleStrike   = probDebrisHere * A_Proj * V_impact * delta_t;
                //                                double probOfCasualty       = probDebrisHere * A_Casualty * V_impact * delta_t;
                //                                double probOfCatastrophe    = probDebrisHere * A_Catastrope * V_impact * delta_t;
                
                double expectedNumPiecesHere        = numberOfPiecesMean[curID];      // This is the number of debris from this debIX in the catalog
                double probOfNoStrikeFromCurID      = pow(1. - probOfSingleStrike, expectedNumPiecesHere);
                //                                double probOfNoCasualtyFromCurID    = pow(1. - probOfCasualty, expectedNumPiecesHere);
                //                                double probOfNoCatastropheFromCurID = pow(1. - probOfCatastrophe, expectedNumPiecesHere);
                
                probNoStrike        *= probOfNoStrikeFromCurID;
                //                                probNoCasualty      *= probOfNoCasualtyFromCurID;
                //                                probNoCatastrophe   *= probOfNoCatastropheFromCurID;
                
            }
            
            // Want to save the time at which the aircraft is first at risk
            if (probabilityOfImpactRecord[acid] == 0.){
                secondsSinceImpact = debris_tx * getDeltaT();
                
            }
            
            double probStrike = (1. - probNoStrike);
            // When initialized, will initialize to 0.
            // Update rule is to multiply the old probNoStrike by the new probNoStrike to get the total probNoStrike
            // Then 1 - that is the probability of a strike
            probabilityOfImpactRecord[acid] = 1- (1-probabilityOfImpactRecord[acid])*probNoStrike;
            
            outfile << "debris_tx = " << debris_tx << ", plane_seconds = " << plane_seconds << ", probStrike = " << probStrike << ", total = " << probabilityOfImpactRecord[acid] << endl;
        }
        cout << "minutes since accident = " << secondsSinceImpact/60. << endl;

    }
    
                        

    
    
    
    //    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
    //
    //        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
    //
    //        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
    //            int zindex = it_z->first;
    //            //cout << "  " << xindex << endl;
    //
    //            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
    //                int xindex = it_x->first;
    //                //cout << "    " << yindex << endl;
    //
    //                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
    //                    int yindex = it_y->first;
    //
    //                    // ========== This is roughly how Wilde computes the probabilities =============
    //                    double probNoStrike         = 1.;
    //                    double probNoCasualty       = 1.;
    //                    double probNoCatastrophe    = 1.;
    
    outfile.close();
    
    return probabilityOfImpactRecord;
}



















map<int, double> SkyGrid::CalculateRiskToIndividualAircraft(vector<int> numberOfPiecesMean, vector<double> arefMean,int secondsFromMidnightUTC){
    // NOTE: I'm going to require that the aircraft data come at a fixed timestep so if you ever wind up
    //   having an irregular timestep then this will break.  Don't break it.
    int planeStepSec = 1;
    
    double delta_t = getDeltaT();
    double ft_2_km = 0.0003048;
    
    string debugFileName = "debugOutput.txt";
    ofstream outfile;
	outfile.open(debugFileName.c_str(), ios::out);
        
    
//    vector<int> pieces (numberOfPiecesMeanArray, numberOfPiecesMeanArray + numDebIX);
//    vector<double> areas (arefMeanList, arefMeanList + numDebIX);
    
    // From python aircraftRecord[acid][0].append([float(curTrackTime), curLat, curLon, curLevel, curSpeed])
    // map<int, pair<vector<vector<double> >, string> >
    // AircraftTrackMap;
    
    // [acType][property] = value
    // property --> dict(topArea=topArea, frontArea=frontArea, acClass=acClass)
    // acClass = 0 for business, 1 for commercial
    //    map<string,map<string,double> > AircraftPropertiesMap;
    
    
    map<int, pair< vector< vector<double> >, string> >::iterator it_ACID;
    
    for (it_ACID = AircraftTrackMap.begin(); it_ACID != AircraftTrackMap.end(); ++it_ACID){
        int acid        = it_ACID->first;           // the aircraft id -- self-explanatory
        string acType   = it_ACID->second.second;   // B737, etc
    
        cout << "acid = " << acid << ", model = " << acType << endl;
        int lenHere = (it_ACID->second.first).size();
        
        // Find the index where the explosion occurs
        int tx = 0;
        while (tx < lenHere){
            if ((it_ACID->second.first)[tx][0] == secondsFromMidnightUTC) {
                break; }
            
            tx++;
        }
        
        for (int ix = 0; (tx + ix) < std::min(lenHere,tx + 5); ix++){
            cout << "   " << (it_ACID->second.first)[tx + ix][0] << "   "
                    << (it_ACID->second.first)[tx + ix][1] << "  "
                    << (it_ACID->second.first)[tx + ix][2] << "  "
                    << (it_ACID->second.first)[tx + ix][3] << "  "
                    << (it_ACID->second.first)[tx + ix][4] << "  "
                    << endl;
        }
    }
    
    map<string,map<string,double> >::iterator it_acType;
    map<string,double>::iterator it_property;
    
    cout << "========= Dumping the AircraftPropertiesMap ========= " << endl;
    for (it_acType = AircraftPropertiesMap.begin(); it_acType != AircraftPropertiesMap.end(); ++it_acType){
        cout << it_acType->first << endl;
        for (it_property = it_acType->second.begin(); it_property != it_acType->second.end(); ++it_property){
            cout << "   " << it_property->first << " = " << it_property->second << endl;
        }
    }
    
    // Start the calculation

    // The iterators for ProbabilityMapDebIX
    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
    map<int, map<int, map<int,binData> > >::iterator it_x;
    map<int, map<int, binData> >::iterator it_y;
    map<int, binData>::iterator it_ID;
    
    double cellVolume = xBinLength*yBinLength*zBinHeight;

    map<int, double> probabilityOfImpactRecord;
    
    // Loop over the aircraft themselves
    // NEED THE TIME OFFSET!!!
    
    for (it_ACID = AircraftTrackMap.begin(); it_ACID != AircraftTrackMap.end(); ++it_ACID){
        int acid        = it_ACID->first;           // the aircraft id -- self-explanatory
        string acType   = it_ACID->second.second;   // B737, etc
        
        probabilityOfImpactRecord[acid] = 0.;
        
        cout << "acid = " << acid << ", model = " << acType << endl;
        int lenHere = (it_ACID->second.first).size();
        
        // Find the areas that are relevant for this aircraft
//        acClass = 1
//        frontArea = 0.0223973
//        topArea = 0.154464
        int     acClass     = AircraftPropertiesMap[acType]["acClass"];
        double  frontArea   = AircraftPropertiesMap[acType]["frontArea"];
        double  topArea     = AircraftPropertiesMap[acType]["topArea"];
        
        // Find the index where the explosion occurs.  Time Explosion Offest = TEO
        int plane_tx = 0;
        while (plane_tx < lenHere){
            if ((it_ACID->second.first)[plane_tx][0] == secondsFromMidnightUTC) {
                break; }
            
            plane_tx++; }
        int TEO = (it_ACID->second.first)[plane_tx][0];
        
        cout << "TEO = " << TEO << endl;
        
        Point tempPt;
        for (plane_tx = plane_tx; plane_tx < lenHere; plane_tx++){
            
            // What's the time index at this point?
            double plane_seconds       = (it_ACID->second.first)[plane_tx][0];
            int debris_tx           = (int) floor((plane_seconds - TEO)/all_points_delta_t);
            double aircraftLat      = (it_ACID->second.first)[plane_tx][1];
            double aircraftLon      = (it_ACID->second.first)[plane_tx][2];
            double aircraftZ        = (it_ACID->second.first)[plane_tx][3];
            double aircraftSpeed    = (it_ACID->second.first)[plane_tx][4];
            
            // This should be the xyz coordinate of the map
            tempPt.set_xy_from_latlon(aircraftLat * PI/180., aircraftLon * PI/180.);
            int xindex = floor((tempPt.get_x() - XREF)/xBinLength);
            int yindex = floor((tempPt.get_y() - YREF)/yBinLength);
            int zindex = floor((aircraftZ - ZREF)/zBinHeight);

            // Does it exist?
//            it_time = ProbabilityMapDebIX.find(debris_tx);
//            it_z    = it_time->second.find(zindex);
//            it_x    = it_z->second.find(xindex);
//            it_y    = it_x->second.find(yindex);
            
            // Get iterators to entry in map that corresponds to the time that this aircraft is in the sky.  It may not exist.
            it_time = ProbabilityMapDebIX.find(debris_tx);
            
            // Does that timestep exist?
            if (it_time != ProbabilityMapDebIX.end()) {
                
                // It does, now look for the iterator to the entry in the map that corresponds to the altitude level specified by zindex
                it_z    = it_time->second.find(zindex);
                
//                cout << it_time->first << endl;
                
                // Does it exist?
                if (it_z != it_time->second.end()) {
                    // It does.  Look for an iterator to the current xvalue.
                    it_x    = it_z->second.find(xindex);
                    
                    // Find the lower and upper bound of x values at this point
                    int lowerX = it_z->second.begin()->first;
                    int upperX = (--(it_z->second.end()))->first;
                    
                    //  Print some stuff so I know I'm not crazy
//                    cout << "debris_tx = " << debris_tx << " = " << it_time->first << ", zindex = " << zindex << " = " << it_z->first << endl;
//                    cout << "   lowerX = " << lowerX << ", xindex = " << xindex << ", upperX = " << upperX << endl;
                
                    
                    if (it_x != it_z->second.end()){
                        
                        // Look for an iterator to the current yvalue
                        it_y    = it_x->second.find(yindex);

                        // Does it exist?
                        if (it_y != it_x->second.end()){
                            
                            // It does.  The plane will probabilistically strike some debris here.
                            // Loop through the possible pieces of debris that might be here and calculate probabilities of impact.
                            
                            double probNoStrike         = 1.;

                            for (it_ID = it_y->second.begin(); it_ID != it_y->second.end(); ++it_ID){
                                int curID = it_ID->first;
                                binData PD = (it_ID->second);
                                // Assuming that if we know there is debris in this cell, the location of that debris is uniformly likely to be anywhere in the volume.
                                
                                double probDebrisHere = (PD.probDebris)/cellVolume;
                                
                                //                        double AF = pow(d_Airplane_front + sqrt(PD.avgArea),2);
                                //                        double AT = pow(d_Airplane_top + sqrt(PD.avgArea),2);
                                //                        double probOfSingleStrike = probDebrisHere * (AF*speed787 + AT*PD.avgVel) * delta_t;
                                
                                // Equation takes mass in grams, outputs ft^2, so need to convert to km^2
                                double theta = atan2(aircraftSpeed, PD.avgVel);  // Sure hope velocity is in km/s
                                double A_Proj = frontArea*sin(theta) + topArea*cos(theta);
                                double V_impact = sqrt(pow(aircraftSpeed, 2) + pow(PD.avgVel,2));

                                double A_Casualty   = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);   // This is the case for pieces over 300g
                                double A_Catastrope = pow( sqrt(A_Proj) + sqrt(PD.avgArea) ,2);

                                if (PD.avgMass < 0.001){
                                    // Throw out pieces that are less than 1g.  They pose no danger
                                    A_Casualty      = 0.;
                                    A_Catastrope    = 0.;
                                }
                                else if (PD.avgMass < 0.300){
                                    A_Casualty      = ((0.0085*pow(PD.avgMass * 1e3,2) + 8.5*(PD.avgMass * 1e3) + 200)) * pow(ft_2_km,2);
                                    A_Catastrope    = (0.025 * pow(PD.avgMass * 1e3,2) + 4*(PD.avgMass * 1e3)) * pow(ft_2_km,2);
                                }
                                
                                double probOfSingleStrike   = probDebrisHere * A_Proj * V_impact * delta_t;
//                                double probOfCasualty       = probDebrisHere * A_Casualty * V_impact * delta_t;
//                                double probOfCatastrophe    = probDebrisHere * A_Catastrope * V_impact * delta_t;
                                
                                double expectedNumPiecesHere        = numberOfPiecesMean[curID];      // This is the number of debris from this debIX in the catalog
                                double probOfNoStrikeFromCurID      = pow(1. - probOfSingleStrike, expectedNumPiecesHere);
//                                double probOfNoCasualtyFromCurID    = pow(1. - probOfCasualty, expectedNumPiecesHere);
//                                double probOfNoCatastropheFromCurID = pow(1. - probOfCatastrophe, expectedNumPiecesHere);
                                
                                probNoStrike        *= probOfNoStrikeFromCurID;
//                                probNoCasualty      *= probOfNoCasualtyFromCurID;
//                                probNoCatastrophe   *= probOfNoCatastropheFromCurID;
                                
                            }

                            double probStrike = (1. - probNoStrike);
                            
                            // When initialized, will initialize to 0.
                            // Update rule is to multiply the old probNoStrike by the new probNoStrike to get the total probNoStrike
                            // Then 1 - that is the probability of a strike
                            probabilityOfImpactRecord[acid] = 1- (1-probabilityOfImpactRecord[acid])*probNoStrike;
                            
                            outfile << "debris_tx = " << debris_tx << ", plane_seconds = " << plane_seconds << ", probStrike = " << probStrike << ", total = " << probabilityOfImpactRecord[acid] << endl;
                        }  
                         
                   
                    }
                    
                }
            }
            
            
////            if (it_time != ProbabilityMapDebIX.end()){
//            if (it_x != it_z->second.end()){
////                cout << "woot!" << endl;
//
//                cout << "   " << it_time->first << endl;
//            }
            
        
            
        }
    }
    


//    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
//
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//
//        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
//            int zindex = it_z->first;
//            //cout << "  " << xindex << endl;
//
//            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//                //cout << "    " << yindex << endl;
//
//                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//
//                    // ========== This is roughly how Wilde computes the probabilities =============
//                    double probNoStrike         = 1.;
//                    double probNoCasualty       = 1.;
//                    double probNoCatastrophe    = 1.;
    
    outfile.close();
    
    return probabilityOfImpactRecord;
}








//// This function is mainly only for debugging purposes
//vector<vector<double> >SkyGrid::SendGridToPython(int tx_desired){
//    
//    vector<vector<double> > GridAsVector;   //[tx][
//
//    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh
//    cout << "NumTimestepsHere = " << NumTimestepsHere << endl;
//
//    // Clear this before we load it up
//    GridAsVector.assign(NumTimestepsHere, vector<double>());
//    
//    int numDebrisCategories = countDebrisIX.size();
//    
//    // First must dump the debIX grids down to a single grid
//    map<int, map<int, map<int, map<int,double> > > >CollapsedMap;
//    
//    map<int, map<int, map<int, map<int, map<int,binData> > > > >::iterator it_time;
//    map<int, map<int, map<int, map<int,binData> > > >::iterator it_z;
//    map<int, map<int, map<int,binData> > >::iterator it_x;
//    map<int, map<int, binData> >::iterator it_y;
//    map<int, binData>::iterator it_ID;
//    
//    // Find the probability vector at this timestep for all points in the grid.
//    //    vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
//    
////    for (it_time=ProbabilityMapDebIX.begin(); it_time != ProbabilityMapDebIX.end(); ++it_time) {
//    
//    it_time = ProbabilityMapDebIX.find(tx_desired);
//    if (it_time != ProbabilityMapDebIX.end()){
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//        
//        
//        // For this area, we need to find the min and max latitudes
//        // These are the indicies which are equally spaced
//        int minx = 1e9;
//        int maxx = -1e9;
//        int miny = 1e9;
//        int maxy = -1e9;
//        for (it_z = ProbabilityMapDebIX[tx].begin(); it_z != ProbabilityMapDebIX[tx].end(); ++it_z){
//            int zindex = it_z->first;
//            
//            // For each z-level, find the first (min) xval
//            it_x = ProbabilityMapDebIX[tx][zindex].begin();
//            minx = std::min(it_x->first, minx);
//            
//            // For each z-level, find the last (max) xval
//            it_x = --(ProbabilityMapDebIX[tx][zindex].end());
//            maxx = std::max(it_x->first, maxx);
//            
//            // For each z-level, find the first (min) yval
//            it_x = ProbabilityMapDebIX[tx][zindex].begin();
//            minx = std::min(it_x->first, minx);
//            
//            // For each z-level, find the last (max) yval
//            it_x = --(ProbabilityMapDebIX[tx][zindex].end());
//            maxx = std::max(it_x->first, maxx);
//        }
//        
//        cout << "minx,max = " << minx << ", " << maxx << endl;
//        
//            
//        
//            
//            
//            
//            
//            
////            for (it_x = ProbabilityMapDebIX[tx][zindex].begin(); it_x != ProbabilityMapDebIX[tx][zindex].end(); ++it_x){
////                int xindex = it_x->first;
////                
////                for (it_y = ProbabilityMapDebIX[tx][zindex][xindex].begin(); it_y != ProbabilityMapDebIX[tx][zindex][xindex].end(); ++it_y){
////                    int yindex = it_y->first;
////                    
////                    double curProb = 0.;
////                    for (it_ID = ProbabilityMapDebIX[tx][zindex][xindex][yindex].begin(); it_ID != ProbabilityMapDebIX[tx][zindex][xindex][yindex].end(); ++it_ID){
////                        int curID = it_ID->first;
////                                                
////                        binData PD = (it_ID->second);
////                        double probDebrisHere = ((it_ID->second).probDebris);
////                        
////                        curProb += probDebrisHere;
////                        
//////                        CollapsedMap[tx][zindex][xindex][yindex] += probDebrisHere/numDebrisCategories;
////                        
////                    }
////                    
////                    // =========== End of the ACTA method ===========================================
////                    
//////                    // Save the probability for this grid cell.  Only put it in the leading debIX.
//////                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probImpact         = probStrike;
//////                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCasualty       = probOfAirplaneInCell*(1. - probNoCasualty)*pFail;
//////                    ProbabilityMapDebIX[tx][zindex][xindex][yindex][firstIX].probCatastrophe    = probOfAirplaneInCell*(1. - probNoCatastrophe)*pFail;
////                    
////                } } }
//    } else {
//        cout << "This time does not exist" << endl;
//    }
//
//    // Then format the grid into a vector that can be passed back
//    
//    
//    
//    
////
//
////    
////    Point tempPt;
//    
//    
//    
//
//
//    
////    // Here is how you iterate through the map
////    map<int, map<int, map<int, map<int,double> > > >::iterator it_time;
////    map<int, map<int, map<int,double> > >::iterator it_z;
////    map<int, map<int,double> >::iterator it_x;
////    map<int,double>::iterator it_y;
////    
////    for (it_time=GridMap.begin(); it_time!=GridMap.end(); it_time++) {
////        
////        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
////        //cout << tx << /*" => " << it_time->second <<*/ '\n';
////        
////        for (it_z = GridMap[tx].begin(); it_z != GridMap[tx].end(); ++it_z){
////            int zindex = it_z->first;
////            //cout << "  " << xindex << endl;
////            
////            for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
////                int xindex = it_x->first;
////                //cout << "    " << yindex << endl;
////                
////                for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
////                    int yindex = it_y->first;
////                    
////                    // So at this point, we're just asking "where are there points?"
////                    //   We don't actually care how many points were present, so we won't even be using the value in GridMap[tx][zindex][xindex][yindex]
////                    //   It's enough that we simply made it here with the iterators
////                    
////                    double lowerleftX = xref + xindex*xBinLength;   // Note that xindex is most likely negative
////                    double lowerleftY = yref + yindex*yBinLength;   // Note that yindex may be positive or negative
////                    double lowerleftZ = zref + zindex*zBinHeight;
////                    
////                    // Want to output the bounds of the non-zero grid cell
////                    for (int dx = 0; dx < 2; dx++){
////                        for (int dy = 0; dy < 2; dy++){
////                            // BIG KLUDGE!!!  Hopefully fixes the fact that the swinging arm doesn't sort left-to-right
////                            double randx = (rand() % 1000)/1000000.;
////                            double randy = (rand() % 1000)/1000000.;
////                            
////                            tempPt.set_xyz(lowerleftX + dx*xBinLength + randx, lowerleftY + dy*yBinLength + randy, lowerleftZ + eps);
////                            
////                            
////                            all_points_total[tx].push_back(tempPt);
////                            
////                        } }
////                    
////                    
////                }
////            }
////            
////            
////            
////        } }
//    
//    
//    
//    return GridAsVector;
//}







//// Times specified here as inputs are in GMT
//// They are sent off to GoogleEarth as GMT but GE MIGHT CHANGE THEM to be in your local time when it renders them (so PST for me)
////  Check whether or not you have that option specified in the GE GUI
//void SkyGrid::ExportBinnedDebrisGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min){
//
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use ExportBinnedDebrisGoogleEarth after converting to probabilities\n\n";
//        exit(-13);
//    }
//
//    // This might be the only function that uses this vector, so just make it here at the end
//    vector<vector<Point> > GridBoundaryPts;
//    GridBoundaryPts.assign(xNumBins+1, vector<Point>() );
//    for (int xx = 0; xx < (xNumBins+1); xx++){
//        GridBoundaryPts[xx].assign(yNumBins+1, Point() );
//        for (int yx = 0; yx < (yNumBins+1); yx++){
//            GridBoundaryPts[xx][yx].set_xyz(xmin + xx * xBinLength, ymin + yx * yBinLength, 0.);
//        }
//    }
//
//
//
//
//
//    // This function assumes a constant timestep delta_t
//    //    void Footprint3D::exportGoogleEarth(char *googleEarthFile, int yyyy, int mm, int dd, int hour, int min){
//    KmlFactory* factory(KmlFactory::GetFactory());
//    kmldom::DocumentPtr document(factory->CreateDocument());
//    document->set_description("This is TJC description right here");
//
//    // SCALING FACTOR!!!!
//    double scaleFactor = 10.;
//
//    // --------------------------------------------------------
//    //Get all the kml stuff set up ----------------------------------------
//    //    KmlFactory* factory(KmlFactory::GetFactory());
//    //    kmldom::DocumentPtr document(factory->CreateDocument());
//    //    document->set_description("This is TJC description right here");
//
//
//    string rocketStyleString = "m_ylw-pushpin";
//    string rocketString = "s_ylw-pushpin";
//    string rocketString_hl = "s_ylw-pushpin_hl";
//    string rocketHref = "http://www.clker.com/cliparts/5/a/8/7/12375609571200265874pitr_Rocket_icon.svg.med.png";
//
//    kmlbase::Color32 lineStyleColor;
//    lineStyleColor.set_color_abgr("ff1e17ff");
//
//    kmlbase::Color32 polyStyleColor;
//    polyStyleColor.set_color_abgr("ff2427ff");
//
//    // Create the StyleMaps -----------------------------------------||
//
//    // Rocket Normal Pair     --------------------
//    kmldom::PairPtr rocketPair( factory->CreatePair());
//    rocketPair->set_key(0);
//    rocketPair->set_styleurl(rocketString);
//
//    // Rocket Highlighted Pair --------------------
//    kmldom::PairPtr rocketPair_hl( factory->CreatePair());
//    rocketPair_hl->set_key(1);
//    rocketPair_hl->set_styleurl(rocketString_hl);
//
//    kmldom::StyleMapPtr rocketStyleMap(factory->CreateStyleMap());
//    rocketStyleMap->set_id(rocketStyleString);
//    rocketStyleMap->add_pair(rocketPair);
//    rocketStyleMap->add_pair(rocketPair_hl);
//
//    document->add_styleselector(rocketStyleMap);
//
//
//
//
//
//
//    // Define the style created in the StyleMap for ROCKET -------------------
//    //    kmldom::IconStyleIconPtr rocketIcon( factory->CreateIconStyleIcon() );
//    //    rocketIcon->set_href(rocketHref);
//    //
//    //    kmldom::IconStylePtr rocketIconStylePtr( factory->CreateIconStyle() );
//    //    rocketIconStylePtr->set_icon(rocketIcon);
//
//    kmldom::LineStylePtr lineStylePtr ( factory->CreateLineStyle() );
//    lineStylePtr->set_color(lineStyleColor);
//
//    kmldom::PolyStylePtr polyStylePtr ( factory->CreatePolyStyle() );
//    polyStylePtr->set_color(polyStyleColor);
//
//
//    kmldom::StylePtr rocketStyle( factory->CreateStyle() );
//    rocketStyle->set_id(rocketString);
//    rocketStyle->set_linestyle(lineStylePtr);
//    rocketStyle->set_polystyle(polyStylePtr);
//    //    rocketStyle->set_iconstyle(rocketIconStylePtr);
//
//    document->add_styleselector(rocketStyle);
//
//
//    // Define the style created in the StyleMap for ROCKET_HL -------------------
//    //    kmldom::IconStyleIconPtr rocketIcon_hl( factory->CreateIconStyleIcon() );
//    //    rocketIcon_hl->set_href(rocketHref);
//    //
//    //    kmldom::IconStylePtr rocketIconStylePtr_hl( factory->CreateIconStyle() );
//    //    rocketIconStylePtr_hl->set_icon(rocketIcon_hl);
//    //
//    //    kmldom::StylePtr rocketStyle_hl( factory->CreateStyle() );
//    //    rocketStyle_hl->set_id(rocketString_hl);
//    //    rocketStyle_hl->set_iconstyle(rocketIconStylePtr_hl);
//    //
//    //    document->add_styleselector(rocketStyle_hl);
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
//    // -------------------------------------------------------
//
//
//
//
//    // Create the initial time structure
//    time_t rawtime;
//    time(&rawtime);
//
//    struct tm * timeinfo;
//    timeinfo = localtime ( &rawtime );  //Have to initialize the structure to the current local time then change it
//    //    timeinfo = gmtime ( &rawtime );  //Have to initialize the structure to the current local time then change it
//
//
//
//    yyyy = yyyy - 1900; //put year into proper format
//    mm = mm - 1;        //how many months SINCE january (so jan mm = 0)
//
//    timeinfo->tm_year = yyyy;
//    timeinfo->tm_mon = mm;
//    timeinfo->tm_mday = dd;
//    timeinfo->tm_hour = hour;
//    timeinfo->tm_min = min;
//    timeinfo->tm_sec = 0;
//
//    char buffer [80];
//    strftime (buffer,80, "%FT%XZ", timeinfo);
//    cout << buffer << endl;
//
//    //turn rawtime into the running local time
//    rawtime = mktime(timeinfo);
//
//    char startTimeBuf[80];
//    char endTimeBuf[80];
//
//    double delta_t = getDeltaT();
//    int numTimeSteps = Grid.size();
//
//    //    // Print it out for debugging
//    //    for (int tx = 0; tx < numTimeSteps; tx++){
//    //        for (int xx = 0; xx < xNumBins; xx++){
//    //            for (int yx = 0; yx < yNumBins; yx++){
//    //                for (int zx = 0; zx < zNumBins; zx++){
//    //                    cout << tx << xx << yx << zx << " = " << Grid[tx][xx][yx][zx] << endl;  } } } }
//
//
//    double KmToMeters = 1e3;
//
//    double alt_base = 0;
//    int counter = 0;
//    // t = timestep, z = z_bin, h = hull, i = a point
//    for (int t = 0; t < numTimeSteps; t++) {
//
//        //Update runningTime as necessary
//        time_t startTime = rawtime + t*delta_t*60.;
//        timeinfo = localtime ( &startTime );
//        strftime (startTimeBuf,80, "%FT%XZ", timeinfo);
//
//        time_t stopTime = rawtime + (t+1)*delta_t*60.;
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
//        double numPtsHere = totalNumPtsAtStep[t];
//        cout << "nuPtsHere = " << numPtsHere << endl;
//        for (int xx = 0; xx < (xNumBins); xx++){
//            for (int yy = 0; yy < (yNumBins); yy++){
//
//                double Lat1 = GridBoundaryPts[xx][yy].get_gdLatDeg();
//                double Lon1 = GridBoundaryPts[xx][yy].get_lonDeg();
//
//                double Lat2 = GridBoundaryPts[xx+1][yy].get_gdLatDeg();
//                double Lon2 = GridBoundaryPts[xx+1][yy].get_lonDeg();
//
//                double Lat3 = GridBoundaryPts[xx+1][yy+1].get_gdLatDeg();
//                double Lon3 = GridBoundaryPts[xx+1][yy+1].get_lonDeg();
//
//                double Lat4 = GridBoundaryPts[xx][yy+1].get_gdLatDeg();
//                double Lon4 = GridBoundaryPts[xx][yy+1].get_lonDeg();
//
//                for (int zz = 0; zz < zNumBins; zz++) {
//                    int cur_alt_lo = (INTxx) floor((alt_base + zz*zBinHeight)*KmToMeters);
//                    int cur_alt_hi = (INTxx) floor((alt_base + zz*zBinHeight + Grid[t][xx][yy][zz]*scaleFactor/numPtsHere)*KmToMeters);
//
//                    // Add in the top (deleted code for side-walls, look in archive if needed)
//                    {
//                        CoordinatesPtr coordinates(factory->CreateCoordinates());
//                        coordinates->add_latlngalt(Lat1, Lon1, cur_alt_hi);       //first point, base alt
//                        coordinates->add_latlngalt(Lat2, Lon2, cur_alt_hi);      //second point, base alt
//                        coordinates->add_latlngalt(Lat3, Lon3, cur_alt_hi);     //second point, top alt
//                        coordinates->add_latlngalt(Lat4, Lon4, cur_alt_hi);    //first point, top alt
//
//                        kmldom::LinearRingPtr linearRing(factory->CreateLinearRing());
//                        linearRing->set_coordinates(coordinates);
//                        kmldom::OuterBoundaryIsPtr outerBound(factory->CreateOuterBoundaryIs());
//                        outerBound->set_linearring(linearRing);
//
//                        // <Polygon>
//                        kmldom::PolygonPtr polygon(factory->CreatePolygon());
//                        polygon->set_outerboundaryis(outerBound);
//                        polygon->set_altitudemode(kmldom::ALTITUDEMODE_RELATIVETOGROUND);
//
//                        multiGeometry->add_geometry(polygon);
//                    }
//                }
//            }
//        }
//
//
//        // I Think this is NOT a memory leak because, even though allocating new memory every time, I'm not losing
//        //track of where it is, I'm simply storing that pointer in the folder a few lines later.
//        PlacemarkPtr placemarkTemp(factory->CreatePlacemark());
//        placemarkTemp->set_geometry(multiGeometry);
//        placemarkTemp->set_timeprimitive(timespan);
//        placemarkTemp->set_styleurl("#" + rocketStyleString);
//
//        document->add_feature(placemarkTemp);
//    }
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
//    kml_file->SerializeToString(&xml_output);   // THIS FUNCTION SUCKS UP A TON OF MEMORY!!!!!
//
//    // Look at the output in the terminal if you want
//    // cout << xml_output << endl;
//
//    // Write the xml to file
//    bool failed = kmlbase::File::WriteStringToFile(xml_output, googleEarthFile);
//    if (failed){
//        cout << "SkyGrid::ExportBinnedDebrisGoogleEarth failed" << endl;
//    }
//
//    return;
//}



















//double SkyGrid::generateAllPointsFromSimpleHistogram(double thresh, int Ntotal, int numEventsSimulated, double pFail = 1){
//
//    // Start by converting to probability
//    double weight = 1.0;
////    ConvertToProbability(weight);
//
//    cout << "Commented out a crucial function which means this is broken.  Exiting." << endl;
//    exit(18);
//
//    double EV_strike = generateAllPointsFromProbability( thresh,  Ntotal,  numEventsSimulated, pFail);
//
//    return EV_strike;
//}





//// Functions for probabilistically handling footprints
//void SkyGrid::GridTheSky(){
//    //total_points_at[tstep][piece of debris]
//
//    // Pick a reference point (x,y,z) = (0,0,0) and i suppose all_points_UTC for now
//    double xref = 0.;
//    double yref = 0.;
//    double zref = 0.;
//
//    //    vector<vector<Point> > total_points_at = getAllPoints();
//
//    // Bin the points
//    int numTimeSteps = all_points_total.size();
//
//    for (int tx = 0; tx < numTimeSteps; tx++){
//        for (int px = 0; px < getNumPointsAtTstep(tx); px++){
//            //            Point curPoint = total_points_at[tx][px];
//            Point curPoint = all_points_total[tx][px];
//
//            int xindex = floor((curPoint.get_x() - xref)/xBinLength);
//            int yindex = floor((curPoint.get_y() - yref)/yBinLength);
//            int zindex = floor((curPoint.get_z() - zref)/zBinHeight);
//
//            // Adjust zindex
//            // If z is negative, place it on the ground, if z is too high, ignore it
//            bool zgood = true;
//            if (curPoint.get_z() > NASkm) {
//                zindex = -1;
//                zgood = false;
//            } else if (curPoint.get_z() < 0) {
//                zindex = 0;
//            }
//
//            // Might also get passed in spurious points (too high to care about), skip them
//            //   only add the points that are currently within the NAS
//            if (zgood) {
//                totalNumPtsAtStepMap[tx] += 1;
//                GridMap[tx][zindex][xindex][yindex] += 1;
//            }
//        }
//    }
//
//
//    //    // Dump some stats
//    //    for (map<int,int>::iterator it=totalNumPtsAtStepMap.begin(); it!=totalNumPtsAtStepMap.end(); it++) {
//    //        std::cout << it->first << " => " << it->second << '\n'; }
//
//    return;
//}


//void SkyGrid::ConvertToProbability(double weight){
//
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use ConvertToProbability after converting to probabilities\nDo Nothing\n";
//        return;
//    }
//
////    int numTimeSteps = Grid.size();
//
//
//    // Here is how you iterate through the map
//    map<int, map<int, map<int, map<int,double> > > >::iterator it_time;
//    map<int, map<int, map<int,double> > >::iterator it_z;
//    map<int, map<int,double> >::iterator it_x;
//    map<int,double>::iterator it_y;
//
//    for (it_time=GridMap.begin(); it_time!=GridMap.end(); ++it_time) {
//        double zTotalSum = 0;
//
//        int numPtsHere = 0;
//        int numDivisions = 0;
//
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//        //cout << tx << /*" => " << it_time->second <<*/ '\n';
//
//        for (it_z = GridMap[tx].begin(); it_z != GridMap[tx].end(); ++it_z){
//            int zindex = it_z->first;
//            //cout << "  " << xindex << endl;
//
//            for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//                //cout << "    " << yindex << endl;
//
//                for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//
//                    double numer = (double) GridMap[tx][zindex][xindex][yindex];    //probably faster to use it_y->second
//                    double denom = (double) totalNumPtsAtStepMap[tx];
//
//                    numPtsHere += numer;
//
//                    ProbabilityMap[tx][zindex][xindex][yindex] = numer/denom;
//                    zTotalSum += numer/denom;
//
//                    numDivisions++;
//
//                } } }
//
//
//        double checkSum = (zTotalSum/weight - 1.);
////        cout << std::scientific;
////        cout << "tx = " << tx << "  numPtsHere = " << numPtsHere << "   totalNum = " << totalNumPtsAtStepMap[tx]
////            << "   totalSum = " << zTotalSum << "    weight = " << weight << "     checkSum = " << checkSum << endl;
////        cout << "           numDivisions = " << numDivisions << "    estimatedMachineEpsilon = " << checkSum/numDivisions << endl;
//        if (fabs(checkSum) > 1e-13) { cout << "YOUR PROBABILITIES ARE MESSED UP!!!!!   " << checkSum << endl;
//            exit(-15);
//        }
//
//    }
//
//
//
//
//
////    // Allocate the space [tx][xx][yx][zx] and compute the values
////    Probability.assign(numTimeSteps, std::vector<std::vector<std::vector<double> > >() );
////    for (int tx = 0; tx < numTimeSteps; tx++){
////        Probability[tx].assign(xNumBins, std::vector<std::vector<double> >());
////        for (int xx = 0; xx < xNumBins; xx++){
////            Probability[tx][xx].assign(yNumBins, std::vector<double>());
////            for (int yx = 0; yx < yNumBins; yx++){
////                Probability[tx][xx][yx].assign(zNumBins, 0.);
////                for (int zx = 0; zx < zNumBins; zx++){
////                    Probability[tx][xx][yx][zx] = weight * ((double) Grid[tx][xx][yx][zx])/((double) totalNumPtsAtStep[tx]);
////                } } } }
//
//    //    // Do some checks
//    //    for (int tx = 0; tx < numTimeSteps; tx++){
//    //        for (int zx = 1; zx < zNumBins; zx++){
//    //            for (int xx = 0; xx < xNumBins; xx++){
//    //                for (int yx = 0; yx < yNumBins; yx++){
//    //                    if (Grid[tx][xx][yx][zx]) {
//    //                        cout << "FLAG = " << Grid[tx][xx][yx][zx] << endl;
//    //                    }
//    //                } } } }
//
//
////    // Do some checks
////    for (int tx = 0; tx < numTimeSteps; tx++){
////        double zTotalSum = 0;
////        for (int zx = 0; zx < zNumBins; zx++){
////            double zSubSum = 0;
////            for (int xx = 0; xx < xNumBins; xx++){
////                for (int yx = 0; yx < yNumBins; yx++){
////                    zSubSum += Probability[tx][xx][yx][zx];
////                } }
////
////            //cout << "zSubSum[" << zx << "] = " << zSubSum << endl;
////            zTotalSum += zSubSum;
////        }
////        //cout << "zTotalSum[" << tx << "] = " << zTotalSum << endl;
////
////        double checkSum = (zTotalSum/weight - 1.);
////        if (fabs(checkSum) > 1e-14) { cout << "YOUR PROBABILITIES ARE MESSED UP!!!!!   " << checkSum << endl; }
////    }
//
//    isProbability = true;
//
//
//    return;
//}




//void SkyGrid::ASH2(int tx, int zindex){
//    // Bivariate Average Shifted Histogram
//
//    // I think the way this should ACTUALLY work is that you pick the h's and then the m's
//    //      derive from that.  This algo does it the other way around.  CHANGE THIS EVENTUALLY.
//
//    int m1 = 5;     // Number of histograms to average over in X direction
//    int m2 = 5;     // Number of histograms to average over in Y direction
//
//    double delta1 = xBinLength;     // Size of initial histogram cells
//    double delta2 = yBinLength;
//
//    // Allocate the weights assuming a bivariate kernel
//    map<int, double> w_m1;          // The actual weights
//    map<int, double> w_m2;
//
//    // Find the normalization constants and then initialize the weight vectors
//    double normalizationConst1 = 0.;
//    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
//        normalizationConst1 += biweightKernel(ix/m1); }
//
//    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
//        w_m1[ix] = m1 * biweightKernel(ix/m1)/normalizationConst1; }
//
//    double normalizationConst2 = 0.;
//    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
//        normalizationConst2 += biweightKernel(ix/m2); }
//
//    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
//        w_m2[ix] = m2 * biweightKernel(ix/m2)/normalizationConst2; }
//
//    // h is the smoothing parametere applied in each direction
//    double h1 = m1 * delta1;
//    double h2 = m2 * delta2;
//
////    map<int, map<int, map<int, map<int,double> > > > ProbabilityMap;
//    map<int, map<int,double> > SubProbabilityMap;
//
////    map<int, map<int, map<int, map<int,double> > > >::iterator it_time;
////    map<int, map<int, map<int,double> > >::iterator it_z;
//
//    map<int, map<int,double> >::iterator it_x;
//    map<int,double>::iterator it_y;
//
//
//    double numHere = 0;
//    // Spread the probabilities around
//    for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
//        int xindex = it_x->first;
//
//        for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
//            int yindex = it_y->first;
//            // I think by virtue of this iterator existing, that means there is non-zero histogram here
//            numHere += GridMap[tx][zindex][xindex][yindex];
//
//            // Loop over all the weights and spread the probabilities around
//            for (int i = 1 - m1; i <= m1 - 1; i++){
//                for (int j = 1 - m2; j <= m2 - 1; j++){
//                    SubProbabilityMap[xindex + i][yindex + j] += GridMap[tx][zindex][xindex][yindex] * w_m1[i] * w_m2[j];
//
//                    //cout << "SubProbabilityMap " << SubProbabilityMap[xindex + i][yindex + j] << endl;
//                }}
//
//
//        } }
//
//    // Renormalize (only loop through the non-zero entries of the new distribution, otherwise you'll collect a lot of zeros)
//    for (it_x = SubProbabilityMap.begin(); it_x != SubProbabilityMap.end(); ++it_x){
//        int xindex = it_x->first;
//        for (it_y = SubProbabilityMap[xindex].begin(); it_y != SubProbabilityMap[xindex].end(); ++it_y){
//            int yindex = it_y->first;
//
//            // Dividing by TOTAL number of points at this tstep (over all z) because I want these probabilities
//            //  to be normalized for at every time step.  This is a little different than the book algorithm
//            SubProbabilityMap[xindex][yindex] =
////                SubProbabilityMap[xindex][yindex] / ( ((double) totalNumPtsAtStepMap[tx]) * h1 * h2);
//                SubProbabilityMap[xindex][yindex] / ( numHere * h1 * h2);
//        } }
//
//
//
//    // I might be done???
//
//
//
//    // Dump stuff to Matlab for debugging
//
//
//    // HARDCORE DEBUGGING ACTION
//
//    map<int, map<int,double> > testProbMap;
//
//    // First, just check the starting histogram
//    double checkSum = 0.;
//    for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
//        int xindex = it_x->first;
//        //cout << "    " << yindex << endl;
//
//        for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
//            int yindex = it_y->first;
//
//            testProbMap[xindex][yindex] = GridMap[tx][zindex][xindex][yindex]/numHere;
//            checkSum += testProbMap[xindex][yindex];
//        }
//    }
//
//    cout << "checkSum is = " << checkSum << endl;
//    cout << "numHere = " << numHere << endl;
//
//
//
//
//
//
//
//
//    //The original grid values
//    cout << "\n\n\nGrid = [\n";
//    for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
//        int xindex = it_x->first;
//
//        for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
//            int yindex = it_y->first;
//
////            cout << "  " << xindex << "  " << yindex << "  " << GridMap[tx][zindex][xindex][yindex] << endl;
//    } }
//    cout << "];\n\n";
//
//    // The new probabilities
//    double probSum = 0.;
//    cout << "\n\n\nProb = [\n";
//    for (it_x = SubProbabilityMap.begin(); it_x != SubProbabilityMap.end(); ++it_x){
//        int xindex = it_x->first;
//        for (it_y = SubProbabilityMap[xindex].begin(); it_y != SubProbabilityMap[xindex].end(); ++it_y){
//            int yindex = it_y->first;
//
////            cout << "  " << xindex << "  " << yindex << "  " << SubProbabilityMap[xindex][yindex] << endl;
//            probSum += SubProbabilityMap[xindex][yindex];
//        } }
//    cout << "];\n\n";
//
//
//
//    double sumW1 = 0.;
//    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
//        sumW1 += w_m1[ix]; }
//    cout << "Do the weights some to m1? = " << sumW1 << endl;
//
//    double sumW2 = 0.;
//    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
//        sumW2 += w_m2[ix]; }
//    cout << "Do the weights some to m2? = " << sumW2 << endl;
//
//    double sumW12 = 0.;
//    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
//        for (int iy = 1 - m2; iy <= m2 - 1; iy++){
//            sumW12 += w_m1[ix] * w_m2[iy]; } }
//    cout << "Do the weights some to m1*m2? = " << sumW12 << endl;
//
//    cout << "Does the new probability sum to one? = " << probSum << endl;
//
//
//    cout << "TODO: Missing adequate buffer region around points?  Actually, i think the std::map might take care of this on its own" << endl;
//
////    // Initialize weights in X direction (w_m1(i))
////    for (int ix = 1 - m1; ix < m1 - 1; ix++){
////
////        double wm1i = biweightKernel(ix/m1)
////
////        double sum
////        for (int jx = 1 - m2; jx < m2 - 1; jx++){
////
////        }
////    }
//
//
//
//
//    return;
//}



//void SkyGrid::ASH2(){
//    // Bivariate Average Shifted Histogram
//    // Based on the book Multivariate Density Estimation by DAVID W. SCOTT
//    // Scott chooses the weights to get a discrete approximation of the density function,
//    //      but I'm choosing them so that I get a mass function.  Could just delete a lot of the normalization stuff actually
//
//    // I think the way this should ACTUALLY work is that you pick the h's and then the m's
//    //      derive from that.  This algo does it the other way around.  CHANGE THIS EVENTUALLY.
//
//
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use ASH2 after converting to probabilities\nDo Nothing\n";
//        exit(-16);
//        return;
//    }
//
//    isProbability = true;
//
//
//    // h is the smoothing parametere applied in each direction
//    double h1 = 20.;    // Assuming that 10km is a natural spreading distance
//    double h2 = 20.;
//
//    double delta1 = xBinLength;     // Size of initial histogram cells
//    double delta2 = yBinLength;
//
//    int m1 = ceil(h1/delta1);     // Number of histograms to average over in X direction
//    int m2 = ceil(h2/delta2);     // Number of histograms to average over in Y direction
//
//    // Update the spreading distance to be consistent
//    h1 = m1 * delta1;
//    h2 = m2 * delta2;
//
////    double delta1 = xBinLength;     // Size of initial histogram cells
////    double delta2 = yBinLength;
////
////    int m1 = 5;     // Number of histograms to average over in X direction
////    int m2 = 5;     // Number of histograms to average over in Y direction
////
////    double h1 = m1 * delta1;
////    double h2 = m2 * delta2;
//
//    // Allocate the weights assuming a bivariate kernel
//    map<int, double> w_m1;          // The actual weights
//    map<int, double> w_m2;
//
//    // Find the normalization constants and then initialize the weight vectors
//    double normalizationConst1 = 0.;
//    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
//        normalizationConst1 += biweightKernel(ix/m1); }
//
//    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
//        w_m1[ix] = m1 * delta1 * biweightKernel(ix/m1)/normalizationConst1; }
//
//    double normalizationConst2 = 0.;
//    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
//        normalizationConst2 += biweightKernel(ix/m2); }
//
//    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
//        w_m2[ix] = m2 * delta2 * biweightKernel(ix/m2)/normalizationConst2; }
//
////    map<int, map<int, map<int, map<int,double> > > > ProbabilityMap;
////    map<int, map<int,double> > SubProbabilityMap;
//
//    map<int, map<int, map<int, map<int,double> > > >::iterator it_time;
//    map<int, map<int, map<int,double> > >::iterator it_z;
//    map<int, map<int,double> >::iterator it_x;
//    map<int,double>::iterator it_y;
//
////    cout << "Start ASH...";
//    // Spread the probabilities around with the weights
//    for (it_time=GridMap.begin(); it_time!=GridMap.end(); it_time++) {
//        int tx = it_time->first;
//
//        double normFactor = ((double) totalNumPtsAtStepMap[tx]) * h1 * h2;
//
//        for (it_z = GridMap[tx].begin(); it_z != GridMap[tx].end(); ++it_z){
//            int zindex = it_z->first;
//
//            for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//
//                for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//                    // By virtue of this iterator existing, that means there is a non-zero histogram here
//
//                    // Loop over all the weights and spread the probabilities around
//                    for (int i = 1 - m1; i <= m1 - 1; i++){
//                        for (int j = 1 - m2; j <= m2 - 1; j++){
//                            ProbabilityMap[tx][zindex][xindex + i][yindex + j] += (GridMap[tx][zindex][xindex][yindex] * w_m1[i] * w_m2[j] / normFactor);
//
//                        }}
//                }}
//        }}
////    cout << "DONE\n";
//
////    // Check the weights
////    double sumW1 = 0.;
////    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
////        sumW1 += w_m1[ix]; }
////    cout << "Do the weights some to m1? = " << sumW1 << endl;
////
////    double sumW2 = 0.;
////    for (int ix = 1 - m2; ix <= m2 - 1; ix++){
////        sumW2 += w_m2[ix]; }
////    cout << "Do the weights some to m2? = " << sumW2 << endl;
////
////    double sumW12 = 0.;
////    for (int ix = 1 - m1; ix <= m1 - 1; ix++){
////        for (int iy = 1 - m2; iy <= m2 - 1; iy++){
////            sumW12 += w_m1[ix] * w_m2[iy]; } }
////    cout << "Do the weights some to m1*m2? = " << sumW12 << endl;
////
////
////    // Check the probabilities of the GridMap
////    for (it_time=GridMap.begin(); it_time!=GridMap.end(); it_time++) {
////        int tx = it_time->first;
////        double checkSum = 0.;
////        int numHere = 0;
////
////        for (it_z = GridMap[tx].begin(); it_z != GridMap[tx].end(); ++it_z){
////            int zindex = it_z->first;
////
////            for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
////                int xindex = it_x->first;
////
////                for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
////                    int yindex = it_y->first;
////
////                    checkSum += (GridMap[tx][zindex][xindex][yindex] / ((double) totalNumPtsAtStepMap[tx])); } } }
////
////        if (fabs(checkSum - 1.) > 1e-13) {
////            cout << "checkSumGRID[" << tx << "] - 1 = " << checkSum - 1. << endl;
////            exit(-15);
////        }
////    }
//
////    cout << "   Start check...";
//    // Check the probabilities (Be safe and keep this in.  Don't be an asshole and delete such an important check)
//    for (it_time=ProbabilityMap.begin(); it_time!=ProbabilityMap.end(); it_time++) {
//        int tx = it_time->first;
//        double checkSum = 0.;
//        int numHere = 0;
//
//        for (it_z = ProbabilityMap[tx].begin(); it_z != ProbabilityMap[tx].end(); ++it_z){
//            int zindex = it_z->first;
//
//            for (it_x = ProbabilityMap[tx][zindex].begin(); it_x != ProbabilityMap[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//
//                for (it_y = ProbabilityMap[tx][zindex][xindex].begin(); it_y != ProbabilityMap[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//
//                    checkSum += (ProbabilityMap[tx][zindex][xindex][yindex]); } } }
//
//        if (fabs(checkSum - 1.) > 1e-10) {
//            cout << "checkSum[" << tx << "] - 1 = " << checkSum - 1. << "    " << checkSum << endl;
//            exit(-15);
//        }
//
//
//    }
////    cout << "DONE\n";
//
//
//
//
////    // Renormalize (only loop through the non-zero entries of the new distribution, otherwise you'll collect a lot of zeros)
////    for (it_x = SubProbabilityMap.begin(); it_x != SubProbabilityMap.end(); ++it_x){
////        int xindex = it_x->first;
////        for (it_y = SubProbabilityMap[xindex].begin(); it_y != SubProbabilityMap[xindex].end(); ++it_y){
////            int yindex = it_y->first;
////
////            // Dividing by TOTAL number of points at this tstep (over all z) because I want these probabilities
////            //  to be normalized for at every time step.  This is a little different than the book algorithm
////            SubProbabilityMap[xindex][yindex] =
////            //                SubProbabilityMap[xindex][yindex] / ( ((double) totalNumPtsAtStepMap[tx]) * h1 * h2);
////            SubProbabilityMap[xindex][yindex] / ( numHere * h1 * h2);
////        } }
//
//
//
//    // I might be done???
//
//
//
//    // Dump stuff to Matlab for debugging
//
//
////    // HARDCORE DEBUGGING ACTION
////
////    map<int, map<int,double> > testProbMap;
////
////    // First, just check the starting histogram
////    double checkSum = 0.;
////    for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
////        int xindex = it_x->first;
////        //cout << "    " << yindex << endl;
////
////        for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
////            int yindex = it_y->first;
////
////            testProbMap[xindex][yindex] = GridMap[tx][zindex][xindex][yindex]/numHere;
////            checkSum += testProbMap[xindex][yindex];
////        }
////    }
////
////    cout << "checkSum is = " << checkSum << endl;
////    cout << "numHere = " << numHere << endl;
////
////    //The original grid values
////    cout << "\n\n\nGrid = [\n";
////    for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
////        int xindex = it_x->first;
////
////        for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
////            int yindex = it_y->first;
////
////            //            cout << "  " << xindex << "  " << yindex << "  " << GridMap[tx][zindex][xindex][yindex] << endl;
////        } }
////    cout << "];\n\n";
////
////    // The new probabilities
////    double probSum = 0.;
////    cout << "\n\n\nProb = [\n";
////    for (it_x = SubProbabilityMap.begin(); it_x != SubProbabilityMap.end(); ++it_x){
////        int xindex = it_x->first;
////        for (it_y = SubProbabilityMap[xindex].begin(); it_y != SubProbabilityMap[xindex].end(); ++it_y){
////            int yindex = it_y->first;
////
////            //            cout << "  " << xindex << "  " << yindex << "  " << SubProbabilityMap[xindex][yindex] << endl;
////            probSum += SubProbabilityMap[xindex][yindex];
////        } }
////    cout << "];\n\n";
////
////
////
////
////    cout << "Does the new probability sum to one? = " << probSum << endl;
//
//
////    cout << "TODO: Missing adequate buffer region around points?  Actually, i think the std::map might take care of this on its own" << endl;
//
//    //    // Initialize weights in X direction (w_m1(i))
//    //    for (int ix = 1 - m1; ix < m1 - 1; ix++){
//    //
//    //        double wm1i = biweightKernel(ix/m1)
//    //
//    //        double sum
//    //        for (int jx = 1 - m2; jx < m2 - 1; jx++){
//    //
//    //        }
//    //    }
//
//
//
//
//    return;
//}







//// At some point need to take into account the curvature of the Earth...the sorting that we do in here
////   relies on the grid cell areas all being equal...so we need a map projection scheme that is equal area.
//double SkyGrid::generateAllPointsFromProbability(double thresh, int Ntotal, int numEventsSimulated, double pFail = 1){
//
//    // Thresh is system-level allowable expected value of plane-debris strike
//    // Ntotal is the total number of timesteps in the system-level simulation
//    // numEventsSimulated is the number of debris-generating events simulated
//    //      This should eventually be stored in SkyGrid when you get a chance
//
//    // I THINK YOU DON'T NEED NTOTAL!!!!  Unioning envelopes will only make you safer.
//
//    double ExpectedValueCollisionSum = 0.;  // EV of collision integrated over all timesteps for the given tfail
//
//
//    double xref = 0.;
//    double yref = 0.;
//    double zref = 0.;
//
//    //cout << "RETURN90\n";
//
//
////    // Start by converting to probability
////    double weight = 1.0;
////    ConvertToProbability(weight);
//
//    // Probability of an airplane (calculated in the KDE function)
//    double probOfAirplane = 0.13989152989;
//
//
//
//    // Nhere
//    int NumTimestepsHere = getNumRange();
//
//
//    // Clear this before we load it up
//    all_points_total.clear();
//    all_points_total.assign(NumTimestepsHere, vector<Point>());
//
//    map<int, map<int, map<int, map<int,double> > > >::iterator it_time;
//    map<int, map<int, map<int,double> > >::iterator it_z;
//    map<int, map<int,double> >::iterator it_x;
//    map<int,double>::iterator it_y;
//
//
//
//
//
//    for (it_time=ProbabilityMap.begin(); it_time != ProbabilityMap.end(); ++it_time) {
//
////        if (it_time == ProbabilityMap.begin()){
////            cout << "RETURN9\n";
////            return 9; }
//
//
//        vector<double> temp4Vec;
//        temp4Vec.assign(4,0.);
//
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//
//        // **Expected Value** Pieces of Debris at each time step (simple average for now)
//        double ExpectedValueNumDebrisHere =  ((double) totalNumPtsAtStepMap[tx]) / ((double) numEventsSimulated);
//
////        vector<Point> tempVec;
//        if (ExpectedValueNumDebrisHere > 0){
//            // Find the probability vector at this timestep for all points in the grid.
//            vector<vector<double> > ProbabilityHere;    // Stores the x index, y index, and probability value for every cell at this tstep and zstep
//
//
//            for (it_z = ProbabilityMap[tx].begin(); it_z != ProbabilityMap[tx].end(); ++it_z){
//                int zindex = it_z->first;
//                //cout << "  " << xindex << endl;
//
//                for (it_x = ProbabilityMap[tx][zindex].begin(); it_x != ProbabilityMap[tx][zindex].end(); ++it_x){
//                    int xindex = it_x->first;
//                    //cout << "    " << yindex << endl;
//
//                    for (it_y = ProbabilityMap[tx][zindex][xindex].begin(); it_y != ProbabilityMap[tx][zindex][xindex].end(); ++it_y){
//                        int yindex = it_y->first;
//
//                        double lowerleftX = xref + xindex*xBinLength;   // Note that xindex is most likely negative
//                        double lowerleftY = yref + yindex*yBinLength;   // Note that yindex may be positive or negative
//                        double lowerleftZ = zref + zindex*zBinHeight;
//
//                        temp4Vec[0] = lowerleftX;
//                        temp4Vec[1] = lowerleftY;
//                        temp4Vec[2] = lowerleftZ;
//                        temp4Vec[3] = ProbabilityMap[tx][zindex][xindex][yindex];
//
//                        ProbabilityHere.push_back(temp4Vec);
//
//                    } } }
//
//
//            // Since the grid cell areas are the same AND assuming uniform aircraft density, sorting the CDF gives same ranking as sorting the EsubC values
//            // using myfunction: compares the third index (the probability) of each vector passed to it
//            // these should be sorted from highest to lowest probability
//            std::sort(ProbabilityHere.begin(), ProbabilityHere.end(), compareXYZ);
//
//            // ~~~~~~~~ Expected Casualty Calculation ~~~~~~~~~~~~~
//            int stopIX = ProbabilityHere.size();      // index in reverse (lowest prob first)
//            double probSum = 0.;            // EsubC = probSum * ProbOfAirplane * ExpectedValueNumDebrisHere
////            double probLevel = thresh / (Ntotal * NumTimestepsHere * probOfAirplane * ExpectedValueNumDebrisHere * pFail);
//            double probLevel = thresh / (NumTimestepsHere * probOfAirplane * ExpectedValueNumDebrisHere * pFail);
//
////            cout << "probLevel = " << probLevel << endl;
////            cout << "numCellsHere = " << stopIX << endl;
////            cout << "thresh = " << thresh << endl;
////            cout << "Ntotal = " << Ntotal << endl;
////            cout << "NumTimestepsHere = " << NumTimestepsHere << endl;
////            cout << "probOfAirplane = " << probOfAirplane << endl;
////            cout << "ExpectedValueNumDebrisHere = " << ExpectedValueNumDebrisHere << endl;
//
//            // Figure out where the last index is
//            while (probSum < probLevel){
//                stopIX--;                   // Subtract 1 before accessing vector for first time
//                probSum += ProbabilityHere[stopIX][3];
//                if (stopIX == 0){ break; }
//                //                cout << "probSum[" << stopIX << "] = " << probSum << "  ";
//            }
//            //            cout << "\n\n\n";
//
//            // Technically, that last index would have put us over the limit.  Add one back.
//            if (stopIX > 0){
//                probSum -= ProbabilityHere[stopIX][3];
//                stopIX++;
//            }
//
//            // EV = P(fail) P(fail at tfail | fail) * P(debrisHere | fail at tfail) * E[numDebrisHere | debrisHere] * P(airplaneHere)
//            //    Note: pFail = P(fail) P(fail at tfail | fail)
//            //          ExpectedValueNumDebrisHere = E[total numDebrisGenerated at tx | fail at tfail]
//            //          E[numDebrisHere | debrisHere] = E[total numDebrisGenerated | fail at tfail] * P(debrisHere | fail at tfail)
//            //                                          = ExpectedValueNumDebrisHere * probSum
//            double ExpectedValueCollision = pFail * probSum * ExpectedValueNumDebrisHere * probOfAirplane;
//            ExpectedValueCollisionSum += ExpectedValueCollision;
//            //cout << "Expected value of collision at this tstep = " << std::scientific << ExpectedValueCollision << endl;
////            cout << "\n";
//
//
//
//            // resize the vector and return it
//            ProbabilityHere.resize(stopIX);
//
//            // Turn it into a points vector and return it
//            double eps = 1e-4;  // Add a little bit to make sure this gets binned properly later
//
//            Point hollowPoint;
//
//
////            tempVec.assign(4*stopIX,Point());
//            for (int ix = 0; ix < stopIX; ix++){
//
////                // z-values are same
////                double curz = zmin + ProbabilityHere[ix][2]*zBinHeight + eps;
//
////                int subIX = 0;
//                for (int dx = 0; dx < 2; dx++){
//                    for (int dy = 0; dy < 2; dy++){
//                        // BIG KLUDGE!!!  Hopefully fixes the fact that the swinging arm doesn't sort left-to-right
//                        double randx = (rand() % 1000)/1000000.;
//                        double randy = (rand() % 1000)/1000000.;
//
//                        double curx = ProbabilityHere[ix][0] + dx*xBinLength + randx;
//                        double cury = ProbabilityHere[ix][1] + dy*yBinLength + randy;
//                        double curz = ProbabilityHere[ix][2] + eps;
//
//                        hollowPoint.set_xyz(curx, cury, curz);
//
////                        int lastIX = ProbabilityMap.end()->first;  // the index of the greatest timestep
////                        int firstIX = std::min(0, ProbabilityMap.begin()->first);  // the index of the least timestep
////                        int firstIXREAL = ProbabilityMap.begin()->first;  // the index of the least timestep
//
//                        //cout << tx << "  " << hollowPoint << "   " << all_points_total.size() << "   " << lastIX << "   "
//                            //<< firstIX << "   " << firstIXREAL << endl;
//                        all_points_total[tx].push_back(hollowPoint);
//
////                        tempVec[4*ix + subIX].set_xyz(curx, cury, curz);
////                        subIX++;
//                    } } }
//
//        }
//        else {
//            cout << "There were no debris points here" << endl;
////            tempVec.assign(0,Point());
//        }
//
//
//
////        // Now store the points (leaving it in this form because i'm anticipating splitting this function up)
////        all_points_total[tx].insert(all_points_total[tx].end(), tempVec.begin(), tempVec.end());
////        tempVec.clear();    // So that we can use it for the next timestep
//
//    }
//
//    //cout << "Total expected value of strike for this tfail envelope = " << ExpectedValueCollisionSum << endl;
//
//    return ExpectedValueCollisionSum;
//}



//void SkyGrid::generateAllPointsFromGrid(){
//
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use generateAllPointsFromGrid after converting to probabilities\n\n";
//        exit(-13);
//    }
//
//    string MatlabDebugFile("GeneratedFiles/MatlabDebugFile.m");
//
//    //    vector< vector<Point> > pointsToWrap;
//    all_points_total.clear();
//
////    int timeSteps = Grid.size();
//
////    if (totalNumPtsAtStepMap.begin()->first != 0){
////        // Will need to handle the cases of >0 and <0 so we iterate through GridMap and put them in the all_points correctly
////        cout << "YOU'VE GOT SOME DEBUGGING TO DO!!!\n\n";
////        cout << totalNumPtsAtStepMap.begin()->first << endl;
////        exit(-14);
////    }
//
//
//    double xref = 0.;
//    double yref = 0.;
//    double zref = 0.;
//
//    // In general, may not have points present at every timestep, so use the max time step.
////    int timeSteps = totalNumPtsAtStepMap.end()->first;  // the index of the greatest timestep
////
//////    vector<Point> tempVec;
////
////    all_points_total.assign(timeSteps, vector<Point>());
//    
//    // Nhere
//    int NumTimestepsHere = getNumRange();   // This is good for all_points_total, but NOT good for scaling the thresh
//    
//    // Clear this before we load it up
//    all_points_total.clear();
//    all_points_total.assign(NumTimestepsHere, vector<Point>());
//
//
//    double eps = 1e-4;  // Add a little bit to make sure this gets binned properly later
//    Point tempPt;
//
//    // Here is how you iterate through the map
//    map<int, map<int, map<int, map<int,double> > > >::iterator it_time;
//    map<int, map<int, map<int,double> > >::iterator it_z;
//    map<int, map<int,double> >::iterator it_x;
//    map<int,double>::iterator it_y;
//
//    for (it_time=GridMap.begin(); it_time!=GridMap.end(); it_time++) {
//
//        int tx = it_time->first;        //Assuming, for the moment, that it starts at tx = 0
//        //cout << tx << /*" => " << it_time->second <<*/ '\n';
//
//        for (it_z = GridMap[tx].begin(); it_z != GridMap[tx].end(); ++it_z){
//            int zindex = it_z->first;
//            //cout << "  " << xindex << endl;
//
//            for (it_x = GridMap[tx][zindex].begin(); it_x != GridMap[tx][zindex].end(); ++it_x){
//                int xindex = it_x->first;
//                //cout << "    " << yindex << endl;
//
//                for (it_y = GridMap[tx][zindex][xindex].begin(); it_y != GridMap[tx][zindex][xindex].end(); ++it_y){
//                    int yindex = it_y->first;
//
//                    // So at this point, we're just asking "where are there points?"
//                    //   We don't actually care how many points were present, so we won't even be using the value in GridMap[tx][zindex][xindex][yindex]
//                    //   It's enough that we simply made it here with the iterators
//
//                    double lowerleftX = xref + xindex*xBinLength;   // Note that xindex is most likely negative
//                    double lowerleftY = yref + yindex*yBinLength;   // Note that yindex may be positive or negative
//                    double lowerleftZ = zref + zindex*zBinHeight;
//
//                    // Want to output the bounds of the non-zero grid cell
//                    for (int dx = 0; dx < 2; dx++){
//                        for (int dy = 0; dy < 2; dy++){
//                            // BIG KLUDGE!!!  Hopefully fixes the fact that the swinging arm doesn't sort left-to-right
//                            double randx = (rand() % 1000)/1000000.;
//                            double randy = (rand() % 1000)/1000000.;
//
//                            tempPt.set_xyz(lowerleftX + dx*xBinLength + randx, lowerleftY + dy*yBinLength + randy, lowerleftZ + eps);
//
//
//                            all_points_total[tx].push_back(tempPt);
//
//                            } }
//
//
//                }
//            }
//
//
//
//        } }
//
//
//
//
////    for (int t = 0; t < timeSteps; t++){
////
////        for (int iz = 0; iz < zNumBins; iz++){
////            tempVec = getGridValues(t, iz);   //Why does this return odd numbers of points?
////
////            // cout << "[t,z,num] = " << t << " " << iz << " " << tempVec.size() << endl;
////            all_points_total[t].insert(all_points_total[t].end(), tempVec.begin(), tempVec.end());
////        }
////    }
//
//    return;
//}




//// ============== THESE FILES WERE NOT REWRITTEN / REPLACED WITH MAPS ======================
//
//// Constructor that loads a saved SkyGrid from a file
//SkyGrid::SkyGrid(string SkyGridFileName)
//: PointCloud()
//{
//    
//    ifstream inFile;
//	inFile.open(SkyGridFileName.c_str(), ios::in | ios::binary);
//    
//    // Write the Timing Info from PointCloud
//    inFile.read((char *) &all_points_UTC, sizeof(all_points_UTC));
//    inFile.read((char *) &all_points_delta_t, sizeof(all_points_delta_t));
//    inFile.read((char *) &all_points_num_range, sizeof(all_points_num_range));
//    
//    // Write the location info from PointCloud
//    inFile.read((char *) &all_points_launchLat, sizeof(all_points_launchLat));
//    inFile.read((char *) &all_points_launchLon, sizeof(all_points_launchLon));
//    inFile.read((char *) &all_points_launchAzimuth, sizeof(all_points_launchAzimuth));
//    
//    // Read in the basic grid size info
//    inFile.read((char *) &xBinLength, sizeof(xBinLength));
//    inFile.read((char *) &yBinLength, sizeof(yBinLength));
//    inFile.read((char *) &zBinHeight, sizeof(zBinHeight));
//    
//    inFile.read((char *) &xNumBins, sizeof(xNumBins));
//    inFile.read((char *) &yNumBins, sizeof(yNumBins));
//    inFile.read((char *) &zNumBins, sizeof(zNumBins));
//    
//    inFile.read((char *) &xmin, sizeof(xmin));
//    inFile.read((char *) &xmax, sizeof(xmax));
//    inFile.read((char *) &ymin, sizeof(ymin));
//    inFile.read((char *) &ymax, sizeof(ymax));
//    inFile.read((char *) &zmin, sizeof(zmin));
//    
//    // Read the number of points at each timestep
//    int numTimeSteps;
//    inFile.read((char *) &numTimeSteps, sizeof(numTimeSteps));
//    totalNumPtsAtStep.assign(numTimeSteps,0);
//    inFile.read((char *) &totalNumPtsAtStep[0], numTimeSteps*sizeof(totalNumPtsAtStep[0]));
//    
//    // Finally, read the actual grid
//    Grid.assign(numTimeSteps, vector< vector< vector<int> > >());
//    for (int t = 0; t < numTimeSteps; t++){
//        Grid[t].assign(xNumBins, vector< vector<int> >());
//        for (int x = 0; x < xNumBins; x++ ){
//            Grid[t][x].assign(yNumBins, vector<int>());
//            for (int y = 0; y < yNumBins; y++){
//                Grid[t][x][y].assign(zNumBins,0);
//                // Not actually sure if this will work, but try it
//                inFile.read((char *) &Grid[t][x][y][0], zNumBins*sizeof(int)); } } }
//    
//    inFile.close();
//    
//    halfBufferCells = ceil(bufferDist/(std::min(xBinLength,yBinLength)));
//    
//    isProbability = false;
//    
//    return;
//}
//
//
//void SkyGrid::StoreGrid(string outFileName){
//    
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use StoreGrid after converting to probabilities\n\n";
//        exit(-13);
//    }
//    
//	ofstream outfile;
//	outfile.open(outFileName.c_str(), ios::out | ios::binary);
//    
//    // Write the Timing Info from PointCloud
//    outfile.write((char *) &all_points_UTC, sizeof(all_points_UTC));
//    outfile.write((char *) &all_points_delta_t, sizeof(all_points_delta_t));
//    outfile.write((char *) &all_points_num_range, sizeof(all_points_num_range));
//    
//    // Write the location info from PointCloud
//    outfile.write((char *) &all_points_launchLat, sizeof(all_points_launchLat));
//    outfile.write((char *) &all_points_launchLon, sizeof(all_points_launchLon));
//    outfile.write((char *) &all_points_launchAzimuth, sizeof(all_points_launchAzimuth));
//    
//    // Write the basic grid size info from SkyGrid
//    outfile.write((char *) &xBinLength, sizeof(xBinLength));
//    outfile.write((char *) &yBinLength, sizeof(yBinLength));
//    outfile.write((char *) &zBinHeight, sizeof(zBinHeight));
//    
//    outfile.write((char *) &xNumBins, sizeof(xNumBins));
//    outfile.write((char *) &yNumBins, sizeof(yNumBins));
//    outfile.write((char *) &zNumBins, sizeof(zNumBins));
//    
//    outfile.write((char *) &xmin, sizeof(xmin));
//    outfile.write((char *) &xmax, sizeof(xmax));
//    outfile.write((char *) &ymin, sizeof(ymin));
//    outfile.write((char *) &ymax, sizeof(ymax));
//    outfile.write((char *) &zmin, sizeof(zmin));
//    
//    // Write the number of points at each timestep
//    int numTimeSteps = totalNumPtsAtStep.size();
//    outfile.write((char *) &numTimeSteps, sizeof(numTimeSteps));
//    outfile.write((char *) &totalNumPtsAtStep[0], numTimeSteps*sizeof(totalNumPtsAtStep[0]));
//    
//    // Finally, write the actual grid
//    for (int t = 0; t < numTimeSteps; t++){
//        for (int x = 0; x < xNumBins; x++ ){
//            for (int y = 0; y < yNumBins; y++){
//                // Not actually sure if this will work, but try it
//                outfile.write((char *) &Grid[t][x][y][0], zNumBins*sizeof(int)); } } }
//    
//    outfile.close();
//    
//    return;
//}
//
//// This is used for ordering the points in SkyGrid::getKDE
//bool myfunction (vector<double> i, vector<double> j) { return (i[2] > j[2]); }
//
//vector<Point> SkyGrid::getKDE(double deltaXY, int tstep, int zstep){
//    
//    // Figure out how many points are necessary
//    int numStepsX = ceil((xmax - xmin)/deltaXY);
//    int numStepsY = ceil((ymax - ymin)/deltaXY);
//    int numPts = numStepsX*numStepsY;
//    
//    xpts_KDE.assign(numPts, 0.);
//    ypts_KDE.assign(numPts, 0.);
//    //    vector<double> ans;
//    //    ans.assign(numPts,0.);
//    
//    // Place the X and Y pts
//    int stepIX = 0;
//    for (int ix = 0; ix < numStepsX; ix++){
//        for (int iy = 0; iy < numStepsY; iy++){
//            xpts_KDE[stepIX] = xmin + ix*deltaXY;
//            ypts_KDE[stepIX] = ymin + iy*deltaXY;
//            //            cout << "x = " << xpts_KDE[stepIX] << "    " << ypts_KDE[stepIX] << endl;
//            stepIX++;
//        }
//    }
//    
//    
//    // Run through once and find all the indices that have at least one point in them
//    list<vector<int> > stuffHere;
//    
//    int iz = zstep;
//    vector<int> tempVec;    // WHY IS THIS VECTOR OF **INT** ?????????
//    tempVec.assign(4,-1);
//    
//    int numDebrisAtThisLevel = 0;
//    for (int ix = 0; ix < xNumBins; ix++){
//        for (int iy = 0; iy < yNumBins; iy++){
//            if (Grid[tstep][ix][iy][iz] > 0){
//                tempVec[0] = xmin + ix*xBinLength;
//                tempVec[1] = ymin + iy*yBinLength;
//                tempVec[2] = iz;
//                tempVec[3] = Grid[tstep][ix][iy][iz];
//                stuffHere.push_back(tempVec);
//                
//                numDebrisAtThisLevel += tempVec[3];
//            } } }
//    
//    if (numDebrisAtThisLevel > 0 ) {
//        
//        // Store as points for later sorting purposes
//        
//        //    // Make the arrays to pass off to gsl multiplication routines
//        //    int numImportantPts = stuffHere.size();
//        //    double *ImportantPts = new double[2*numImportantPts];
//        //
//        //    double *xVec = new double[numImportantPts];
//        //    double *yVec = new double[numImportantPts];
//        //
//        //    for (list<vector<int> >::iterator it = stuffHere.begin(); it != stuffHere.end(); it++){
//        //        static int ix = 0;
//        //        ImportantPts[ix] = ((vector<int>) *it)[0];
//        //        ix++;
//        //        ImportantPts[ix] = ((vector<int>) *it)[1];
//        //        ix++;
//        //    }
//        //
//        //    gsl_matrix_view A = gsl_matrix_view_array(ImportantPts, numImportantPts, 2);
//        
//        // Not using a real bandwidth matrix, but pretend
//        double bandwidth = .5;
//        
//        vector< vector<double> > bigOlVector;
//        bigOlVector.assign(numPts, vector<double>() );
//        
//        double cumSum = 0;
//        for (int ix = 0; ix < numPts; ix++){
//            bigOlVector[ix].assign(3,0.);
//            bigOlVector[ix][0] = xpts_KDE[ix];
//            bigOlVector[ix][1] = ypts_KDE[ix];
//            // Place the X and Y pts
//            
//            double expSum = 0;
//            for (list<vector<int> >::iterator it = stuffHere.begin(); it != stuffHere.end(); it++){
//                // (x - xhist)^2 + (y - yhist)^2
//                double innerProd = pow((xpts_KDE[ix] - ((vector<int>) *it)[0]),2) + pow((ypts_KDE[ix] - ((vector<int>) *it)[1]),2);
//                expSum += (((vector<int>) *it)[3]) * exp(-0.5*innerProd/bandwidth);
//            }
//            
//            double answer = expSum/(2*PI*numDebrisAtThisLevel*bandwidth);
//            bigOlVector[ix][2] = answer;
//            cumSum += answer*deltaXY*deltaXY;  //checking that everything integrates to 1
//        }
//        
//        //    // Debug with Matlab
//        //    cout << "\n\n\n\n\n\n\n\n" << endl;
//        //    cout << "close all; clear all; clc;\n\n";
//        //    cout << "xykde = [\n";
//        //    for (int ix = 0; ix < (numPts-1); ix++){
//        //        cout << bigOlVector[ix][0] << "   " << bigOlVector[ix][1] << "   " << bigOlVector[ix][2] << ";\n"; }
//        //    cout << bigOlVector[(numPts-1)][0] << "   " << bigOlVector[(numPts-1)][1] << "   " << bigOlVector[(numPts-1)][2] << "];\n\n\n";
//        
//        
//        // Since the XY step sizes are the same, sorting the PDF gives same ranking as sorting the CDF
//        // using function as comp
//        std::sort(bigOlVector.begin(), bigOlVector.end(), myfunction);
//        
//        
//        
//        // For a 787
//        // typical cruising speed is 900km/h
//        // speed = 954km/h = 15.9km/min
//        // span = 197.31ft  =  0.0601401km
//        // length = 185.56ft = 0.0565587km
//        // top area to avoid = span * (speed*min + length) = 0.95962903587 km^2
//        // 2 n.m. = 3.704 km
//        // 2 (n.m.)^2 = 6.859808 km^2
//        // assume one airplane every 2 n.m.
//        //  ==> probOfAirplane = areaToAvoid/assumedAreaDensity = 0.95962903587/6.859808 = 0.13989152989
//        // 4 n.m. = 7.408 km
//        // 4 (n.m.)^2 * (1.852 km/nm)^2 = 13.72 km^2
//        // assume one airplane ever 4 n.m.
//        //  ==> probOfAirplane = areaToAvoid/assumedAreaDensity = 0.95962903587/13.72 = 0.06994
//
//        double probOfAirplane = 0.13989152989;
//        
//        //    int stopIX = -1;    //prevents overcounting
//        //    double probOfHit = 10;
//        //    double probLevel = 1e-6;
//        //
//        //    cout << "processed = [\n";
//        //    while (probOfHit > probLevel){
//        //        stopIX++;
//        //        probOfHit = bigOlVector[stopIX][2]*deltaXY*deltaXY * probOfAirplane;
//        //        cout << bigOlVector[stopIX][0] << "   " << bigOlVector[stopIX][1] << "   " << bigOlVector[stopIX][2] << "    " << probOfHit << ";\n"; }
//        //    cout << "];\n\n\n";
//        
//        // ~~~~~~~~ Expected Casualty Calculation ~~~~~~~~~~~~~
//        int stopIX = numPts;    // index in reverse
//        double EsubC = 0.;
//        double probLevel = 1e-7;
//        double twoNMsqr = 6.859808; // 2 (nautical miles)^2 in km^2
//        
//        // Figure out where the last index is
//        while (EsubC < probLevel){
//            stopIX--;
//            EsubC += bigOlVector[stopIX][2]*deltaXY*deltaXY * probOfAirplane / (twoNMsqr*twoNMsqr*zNumBins);}   //divide by zNumBins because that many altitude levels
//        
//        //    // output the remaining points
//        //    cout << "processed = [\n";
//        //    for (int ix = 0; ix < stopIX; ix++){
//        //        EsubC += bigOlVector[ix][2]*deltaXY*deltaXY * probOfAirplane / (twoNMsqr*twoNMsqr);
//        //        cout << bigOlVector[ix][0] << "   " << bigOlVector[ix][1] << "   " << bigOlVector[ix][2] << "    " << EsubC << ";\n"; }
//        //    cout << "];\n\n\n";
//        //    // ~~~~~~~~~ End of EsubC Calculation ~~~~~~~~~~~~~~~~~
//        //
//        //
//        //    cout << "numX = " << numStepsX << endl;
//        //    cout << "numY = " << numStepsY << endl;
//        //    cout << "xmesh = reshape(xykde(:,1),numY, numX);\nymesh = reshape(xykde(:,2),numY, numX);\nzmesh = reshape(xykde(:,3),numY, numX);\n\n";
//        //    cout << "contour(xmesh, ymesh, zmesh)\nhold on\nscatter(processed(:,1), processed(:,2), processed(:,4))\n\n";
//        
//        cout << "The probabilities should sum to one = " << cumSum << endl;
//        
//        if (isnan(cumSum)){
//            cout << "hopping into nan debug" << endl;
//        }
//        
//        // resize the vector and return it
//        bigOlVector.resize(stopIX);
//        
//        // Turn it into a points vector and return it
//        double eps = 1e-4;  // Add a little bit to make sure this gets binned properly later
//        
//        vector<Point> ans;
//        ans.assign(stopIX,Point());
//        for (int ix = 0; ix < stopIX; ix++){
//            // BIG KLUDGE!!!  Hopefully fixes the fact that the swinging arm doesn't sort left-to-right
//            double randx = (rand() % 1000)/1000000.;
//            double randy = (rand() % 1000)/1000000.;
//            ans[ix].set_xyz(bigOlVector[ix][0] + randx, bigOlVector[ix][1] + randy, zstep * zBinHeight + eps); }
//        
//        return ans;
//        
//        
//    } else {
//        cout << "There were no debris points here" << endl;
//        
//        vector<Point> ans;
//        ans.assign(0,Point());
//        return ans;
//    }
//}
//
//void SkyGrid::generateAllPointsFromKDE(double deltaXY){
//    
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use generateAllPointsFromKDE after converting to probabilities\n\n";
//        exit(-13);
//    }
//    
//    //    vector< vector<Point> > pointsToWrap;
//    all_points_total.clear();
//    
//    int timeSteps = Grid.size();
//    //    int timeSteps = 4;
//    
//    vector<Point> tempVec;
//    all_points_total.assign(timeSteps, vector<Point>());
//    for (int t = 0; t < timeSteps; t++){
//        for (int iz = 0; iz < zNumBins; iz++){
//            tempVec = getKDE(deltaXY, t, iz);
//            all_points_total[t].insert(all_points_total[t].end(), tempVec.begin(), tempVec.end());
//        }
//    }
//    
//    return;
//}
//
//
//
//
//
//
//
//
//void SkyGrid::ConvertToEmptyProbability(){
//    
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use ConvertToProbability after converting to probabilities\nDo Nothing\n";
//        return;
//    }
//    
//    int numTimeSteps = Grid.size();
//    
//    // Allocate the space [tx][xx][yx][zx] and compute the values
//    Probability.assign(numTimeSteps, std::vector<std::vector<std::vector<double> > >() );
//    for (int tx = 0; tx < numTimeSteps; tx++){
//        Probability[tx].assign(xNumBins, std::vector<std::vector<double> >());
//        for (int xx = 0; xx < xNumBins; xx++){
//            Probability[tx][xx].assign(yNumBins, std::vector<double>());
//            for (int yx = 0; yx < yNumBins; yx++){
//                Probability[tx][xx][yx].assign(zNumBins, 0.);    } } }
//    
//    isProbability = true;
//    
//    return;
//}
//
//
//
//// So we know that in order to do the weighted combinging thing, we will have to already have gotten all the points in their
////  respective skygrids.  Before converting to probabilities, combine all of the skygrids into a single grid like i've been
////  working on for the last few days and use the information from that grid to let us know how many total timesteps and
////  what the grid boundaries are.  This will simplify a lot of the logic and i can just create the correct size of probability matrx
////  right off the bat.
//
//void SkyGrid::weightedCombine(SkyGrid *in, double weight){
//    
//    // Start by making sure everything is a probability
//    if (!isProbability){
//        cout << "weightedCombine ERROR:  You must start by converting the grid to a probability\n\n";
//        exit(-14);
//    }
//    
//    // Make sure certain aspects of both probability grids are the same numbins, sizes, that they nest well, etc
//    cout << "~~~~~~~~~~~~ You owe this function some error checking\n";
//    bool binSizesGood = (xBinLength == in->xBinLength) && (yBinLength == in->yBinLength) && (zBinHeight == in->zBinHeight);
//    bool xBoundsGood = (xmin <= in->xmin) && (xmax >= in->xmax);
//    bool yBoundsGood = (ymin <= in->ymin) && (ymax >= in->ymax);
//    bool xSnapped = fabs(fmod(xmin, xBinLength) - fmod(in->xmin, xBinLength)) <= 1e-14;
//    bool ySnapped = fabs(fmod(ymin, yBinLength) - fmod(in->ymin, yBinLength)) <= 1e-14;
//    
//    cout << "bools are " << binSizesGood << xBoundsGood << yBoundsGood << xSnapped << ySnapped << endl;
//    
//    // Convert the incoming with the desired weight.  If it's already been done, that's fine, nothing will happen.
//    in->ConvertToProbability(weight);
//    
//    // Find the offsets in x, y, and time from the underlying grid (SHOULD ALREADY BE LARGE ENOUGH)
//    int xOffset = round((in->xmin - xmin)/xBinLength);
//    int yOffset = round((in->ymin - ymin)/yBinLength);
//    
//    double curTstepMinutes = getDeltaT();
//    int tOffset = round((in->getInitialUTC() - getInitialUTC())*24*60/curTstepMinutes);
//    
//    // Copy the incoming into the old
//    int inTimeSteps = in->Probability.size();
//    int inXNumBins = in->xNumBins;
//    int inYNumBins = in->yNumBins;
//    
//    for (int tx = 0; tx < inTimeSteps; tx++){
//        for (int xx = 0; xx < inXNumBins; xx++){
//            for (int yy = 0; yy < inYNumBins; yy++){
//                for (int zz = 0; zz < zNumBins; zz++ ){
//                    // Already applied the weight earlier, so no need for that here
//                    Probability[tx + tOffset][xx + xOffset][yy + yOffset][zz] += in->Probability[tx][xx][yy][zz];
//                }
//            } } }
//    
//    return;
//}
//
//
//bool SkyGrid::isTotalProbabilityGood(){
//    
//    // Start by making sure everything is a probability
//    if (!isProbability){
//        cout << "CheckTotalProbability ERROR:  You must start by converting the grid to a probability\n\n";
//        exit(-14);
//    }
//    
//    bool isGood = true;
//    double TotalSum = 0;
//    
//    int numTimeSteps = Probability.size();
//    for (int tx = 0; tx < numTimeSteps; tx++){
//        double zSum = 0;
//        for (int zx = 0; zx < zNumBins; zx++){
//            double zSubSum = 0;
//            for (int xx = 0; xx < xNumBins; xx++){
//                for (int yx = 0; yx < yNumBins; yx++){
//                    zSubSum += Probability[tx][xx][yx][zx];
//                } }
//            
//            //cout << "zSubSum[" << zx << "] = " << zSubSum << endl;
//            zSum += zSubSum;
//        }
//        cout << "zTotalSum[" << tx << "] = " << zSum << endl;
//        TotalSum += zSum;
//    }
//    
//    double checkSum = (TotalSum - 1.);
//    if (fabs(checkSum) > 1e-14) {
//        cout << "YOUR PROBABILITIES ARE MESSED UP!!!!!   " << "  " << TotalSum << endl;
//        isGood = false;
//    }
//    
//    return isGood;
//}
//
//
//// This function feels really sloppy.  Hope to delete it eventually.
//void SkyGrid::RemoveHistogramKeepTiming(){
//    
//    // Check that you can still use it
//    if (isProbability){
//        cout << "ERROR, youre trying to use RemoveHistogramKeepTiming after converting to probabilities\n\n";
//        exit(-13);
//    }
//    
//    int numTimeSteps = totalNumPtsAtStep.size();
//    totalNumPtsAtStep.clear();
//    totalNumPtsAtStep.assign(numTimeSteps, 0);
//    
//    // Finally, read the actual grid
//    for (int t = 0; t < numTimeSteps; t++){
//        for (int x = 0; x < xNumBins; x++ ){
//            for (int y = 0; y < yNumBins; y++){
//                for (int z = 0; z < zNumBins; z++){
//                    Grid[t][x][y][z] = 0;
//                }
//                //                //overwriting
//                //                Grid[t][x][y].clear();
//                //                Grid[t][x][y].assign(zNumBins,0);
//            } } }
//}


















//    // Print it out for debugging
//    for (int tx = 1; tx < numTimeSteps; tx++){
//        for (int xx = 0; xx < xNumBins; xx++){
//            for (int yx = 0; yx < yNumBins; yx++){
//                for (int zx = 0; zx < zNumBins; zx++){
//                    cout << tx << xx << yx << zx << " = " << Grid[tx][xx][yx][zx] << endl;  } } } }
