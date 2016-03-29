//
//  Grid3D.cpp
//  SUFARM
//
//  Created by Thomas Colvin on 3/28/16.
//  Copyright (c) 2016 Thomas Colvin. All rights reserved.
//

#include "Grid3D.h"


Grid3D::Grid3D(){
    SpatialProbabilty[0][0][0] = 666.;
}

/*! Takes a spatial probability grid, stores it
 *
 !*/
Grid3D::Grid3D(map<int, map<int, map<int,double> > > SpatialProbabilty_in){
    // TODO: Eventually can just load the values individially into Grid3D instead of using SpatialProbability_in,
    //   however for now, fastest path forward is just copying
    SpatialProbabilty = SpatialProbabilty_in;
}

Grid3D::Grid3D( const Grid3D &obj){
    SpatialProbabilty = obj.SpatialProbabilty;
}// copy constructor


map<int, map<int, map<int,double> > > Grid3D::getGrid(){
    return SpatialProbabilty;
}
