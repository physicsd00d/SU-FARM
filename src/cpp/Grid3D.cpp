//
//  Grid3D.cpp
//  SUFARM
//
//  Created by Thomas Colvin on 3/28/16.
//  Copyright (c) 2016 Thomas Colvin. All rights reserved.
//

#include "Grid3D.h"


Grid3D::Grid3D(){
//    SpatialProbabilty[0][0][0] = 666.;
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
    printf("Copy Constructor!\n");
    SpatialProbabilty = obj.SpatialProbabilty;
}// copy constructor


map<int, map<int, map<int,double> > > Grid3D::getGrid(){
    return SpatialProbabilty;
}

Grid3D Grid3D::operator+(const Grid3D &in){
//    printf("Plus Operator!\n");
    Grid3D ans;
    ans = in.SpatialProbabilty;      // Start with the incoming grid and add what we've already got to it
    
    map<int, map<int, map<int,double> > >::iterator zit;
    map<int, map<int,double> >::iterator xit;
    map<int,double>::iterator yit;
    
    for (zit = this->SpatialProbabilty.begin(); zit != this->SpatialProbabilty.end(); ++zit){
        int zindex = zit->first;
        
        for (xit = zit->second.begin(); xit != zit->second.end(); ++xit){
            int xindex = xit->first;
            
            for (yit = xit->second.begin(); yit != xit->second.end(); ++yit){
                int yindex = yit->first;
                
                ans.SpatialProbabilty[zindex][xindex][yindex] += yit->second;
            }
        }
    }
    
    return ans;
}

Grid3D Grid3D::operator*(double k){
//    printf("Scalar Mult Operator!\n");
    map<int, map<int, map<int,double> > > ans;
//    ans = in.SpatialProbabilty;      // Start with the incoming grid and add what we've already got to it
    
    map<int, map<int, map<int,double> > >::iterator zit;
    map<int, map<int,double> >::iterator xit;
    map<int,double>::iterator yit;
    
    for (zit = this->SpatialProbabilty.begin(); zit != this->SpatialProbabilty.end(); ++zit){
        int zindex = zit->first;
        
        for (xit = zit->second.begin(); xit != zit->second.end(); ++xit){
            int xindex = xit->first;
            
            for (yit = xit->second.begin(); yit != xit->second.end(); ++yit){
                int yindex = yit->first;
                
                ans[zindex][xindex][yindex] = yit->second * k;
            }
        }
    }
    
    return Grid3D(ans);
}



