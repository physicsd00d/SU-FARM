//
//  Grid3D.h
//  SUFARM
//
//  Created by Thomas Colvin on 3/28/16.
//  Copyright (c) 2016 Thomas Colvin. All rights reserved.
//

#ifndef __SUFARM__Grid3D__
#define __SUFARM__Grid3D__

#include <iostream>

#include <map> // header file needed for to use MAP STL
using std::map;

class Grid3D{
	
private:
    // These would be useful for some error checking, but aren't critical to start
//    double xBinLength, yBinLength, zBinHeight;  //[km]
//    double all_points_UTC;
//    whichProb
    map<int, map<int, map<int,double> > > SpatialProbabilty;

	
public:
    Grid3D();
    Grid3D( const Grid3D &obj);  // copy constructor
    Grid3D(map<int, map<int, map<int,double> > > SpatialProbabilty_in);
	
    map<int, map<int, map<int,double> > > getGrid();
    
    // Overloading operators
    Grid3D operator+(const Grid3D &obj);
    Grid3D operator*(double k);
    bool operator<=(const Grid3D &obj);


};






#endif /* defined(__SUFARM__Grid3D__) */
