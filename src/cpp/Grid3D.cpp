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
    SpatialProbabilty = obj.SpatialProbabilty;
}// copy constructor


map<int, map<int, map<int,double> > > Grid3D::getGrid() const{
    return SpatialProbabilty;
}

Grid3D Grid3D::operator+(const Grid3D &in){
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
    map<int, map<int, map<int,double> > > ans;
    
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

bool Grid3D::operator<=(const Grid3D &obj){
    // this <= obj
    // thus should iterate over this, because if an [z,x,y] isn't in this, then it's automatically 0 <= obj,
    //   assuming obj >= 0 which should be true.  Perhaps should check that as well.
    
    map<int, map<int, map<int,double> > >::iterator zit;
    map<int, map<int,double> >::iterator xit;
    map<int,double>::iterator yit;
    
    for (zit = this->SpatialProbabilty.begin(); zit != this->SpatialProbabilty.end(); ++zit){
        int zindex = zit->first;
        
        for (xit = zit->second.begin(); xit != zit->second.end(); ++xit){
            int xindex = xit->first;
            
            for (yit = xit->second.begin(); yit != xit->second.end(); ++yit){
                int yindex = yit->first;
                if (yit->second > obj.SpatialProbabilty.at(zindex).at(xindex).at(yindex)) {
                    return false;
                }
            }
        }
    }
    
    return true;
}


//Grid3D Grid3D::removeNoDanger(const Grid3D &obj){
//    // For all of the elements in *this, check to see if they have a value in obj.
//    // If they do, then don't touch that element
//    // If they don't have a value, that means they're zero, so set them to zero in *this as well.
//    int ix = 0;
//    map<int, map<int, map<int,double> > > ans = this->getGrid();
//
//    map<int, map<int, map<int,double> > >::iterator zit;
//    map<int, map<int,double> >::iterator xit;
//    map<int,double>::iterator yit;
//    
//    for (zit = ans.begin(); zit != ans.end(); ++zit){
//        int zindex = zit->first;
//        
//        // Make sure there are values at this z-level
//        if (obj.SpatialProbabilty.count(zindex) > 0) {
//        
//            for (xit = zit->second.begin(); xit != zit->second.end(); ++xit){
//                int xindex = xit->first;
//                
//                // Make sure there are values at this x-level
//                if (obj.SpatialProbabilty.at(zindex).count(xindex) > 0) {
//                    
//                    for (yit = xit->second.begin(); yit != xit->second.end(); ++yit){
//                        int yindex = yit->first;
//                        
//                        int count = (int) obj.SpatialProbabilty.at(zindex).at(xindex).count(yindex);
//                        if (count == 0) {
//                            // There is no danger at this location in obj, so set the danger in *this to zero
//                            printf("   [%d][%d][%d] is not present\n", zindex, xindex, yindex);
//                            xit->second.erase(yindex);  // Remove the key completely
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//    return Grid3D(ans);
//    
//}


Grid3D Grid3D::removeNoDanger(const Grid3D &obj){
    // For all of the elements in *this, check to see if they have a value in obj.
    // If they do, then don't touch that element
    // If they don't have a value, that means they're zero, so set them to zero in *this as well.
    map<int, map<int, map<int,double> > > ans = this->getGrid();
    
    map<int, map<int, map<int,double> > >::iterator zit;
    map<int, map<int,double> >::iterator xit;
    map<int,double>::iterator yit;
    
    for (zit = ans.begin(); zit != ans.end(); ++zit){
        int zindex = zit->first;

        // Make sure there are values at this z-level
        int countz = (int) obj.SpatialProbabilty.count(zindex);
        
            for (xit = zit->second.begin(); xit != zit->second.end(); ++xit){
                int xindex = xit->first;
                    
                    for (yit = xit->second.begin(); yit != xit->second.end(); ++yit){
                        int yindex = yit->first;
                        
                        // If there are points at z, check if there are points at x too
                        int countx = 0;
                        if (countz > 0) { countx = (int) obj.SpatialProbabilty.at(zindex).count(xindex); }
                        
                        // There are points at x, check if there are points at y too
                        int county = 0;
                        if (countx > 0){ county = (int) obj.SpatialProbabilty.at(zindex).at(xindex).count(yindex); }
                        
                        if (county == 0) {
                            // There is no danger at this location in obj, so set the danger in *this to zero
                            //printf("   [%d][%d][%d] is not present\n", zindex, xindex, yindex);
                            xit->second.erase(yindex);  // Remove the key completely
                        }
                    }
            }
    }
    
    return Grid3D(ans);
    
}























