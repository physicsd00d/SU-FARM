/*
 *  Point.cpp
 *  Traj_Prop
 *
 *  Created by Thomas Colvin on 8/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include "Point.h"

int Point::sortBy = -1; //Initialize the static variable


ostream& operator<<(ostream& output, const Point& p){
	output << "(" << p.x << ", " << p.y << ", " << p.z << ", " << p.R_local << ")\n";
    output << "  (" << p.Vx_km_s << ", " << p.Vy_km_s << ", " << p.Vz_km_s << ", " << p.mass_kg << ", " << p.area_km2 << ")\n";
	return output;
}



bool Point::operator== (Point rhs) {
	bool isX = (x == rhs.get_x());
	bool isY = (y == rhs.get_y());
	bool isZ = (z == rhs.get_z());
	// Maybe check for local radius too?
	
	return (isX && isY && isZ);
}

Point::Point(){
    x = 0;
    y = 0;
    z = 0;
    R_local = -5.;
    
    Vx_km_s = -99;
    Vy_km_s = -99;
    Vz_km_s = -99;
    mass_kg = -99;
    area_km2 = -99;
    
    calcRefRadius();
}

Point Point::operator- (Point rhs) {
	Point temp;
	
	temp.set_x(x - rhs.get_x());
	temp.set_y(y - rhs.get_y());
	temp.set_z(z - rhs.get_z());

	return temp;
}

Point Point::operator+ (Point rhs) {
	Point temp;
	
	temp.set_x(x + rhs.get_x());
	temp.set_y(y + rhs.get_y());
	temp.set_z(z + rhs.get_z());
	
	return temp;
}

bool Point::operator < (Point rhs) {
    if (sortBy == 1){
        return (x < rhs.get_x()); }
    else if (sortBy == 2){
        return (y < rhs.get_y()); }
    else if (sortBy == 3){
        return (z < rhs.get_z()); }
    else {
        cout << "GOTTA SET THE SORTBY OPTION!!!!~~~~~~~~~~~~~~~~~~~~~~~~~~";
        return false; }
}

double Point::calcAzimuth(const Point &ptFinal){
    // this ptr should point to the initial point we draw the ray from
    // ptFinal should be the point the ray is drawn towards
    // This function finds the azimuth of the ray so defined
    
    double dx = ptFinal.x - this->get_x();
    double dy = ptFinal.y - this->get_y();
    double azimuth = (PI/2. - atan2(dy, dx));
    
    return azimuth;
}

void Point::rotateThisAboutAnotherPoint(const Point &refPoint, double rotAngleRad){
    
    double relX = x - refPoint.x;
    double relY = y - refPoint.y;
    
    double delX = relX*cos(rotAngleRad) - relY*sin(rotAngleRad);
    double delY = relX*sin(rotAngleRad) + relY*cos(rotAngleRad);
    
    x = refPoint.x + delX;
    y = refPoint.y + delY;
    
    return;
}

//Point::Point(double x_in, double y_in, double z_in, double R_local_in) {
//	x = x_in;
//	y = y_in;
//	z = z_in;
//	R_local = R_local_in;
//	num_id = -99;
//}



Point::Point(double gdLatIN, double lonIN, double z_in, double R_local_in) {
	z = z_in;
	R_local = R_local_in;
	num_id = -99;
    //    x = R_equator * lonIN;  // R_equator * longitude_angle
    //    y = R_local * gdLatIN;                // Local Radius of Earth * gdlat

    Vx_km_s = -99;
    Vy_km_s = -99;
    Vz_km_s = -99;
    mass_kg = -99;
    area_km2 = -99;

    calcRefRadius();
//    double a = R_polar * cos(refLat);
//    double b = R_equator * sin(refLat);
//    refRadius = R_equator * R_polar / sqrt(a*a + b*b);
    
    set_xy_from_latlon(gdLatIN, lonIN);
}


Point::Point(const Point &object) {
	x = object.x;
	y = object.y;
	z = object.z;
	R_local = object.R_local;
	num_id = object.num_id;
    
    Vx_km_s = object.Vx_km_s;
    Vy_km_s = object.Vy_km_s;
    Vz_km_s = object.Vz_km_s;
    mass_kg = object.mass_kg;
    area_km2 = object.area_km2;
    
    refRadius = object.refRadius;
}

void Point::calcRefRadius(){
    double a = R_polar * cos(refLat);
    double b = R_equator * sin(refLat);
    refRadius = R_equator * R_polar / sqrt(a*a + b*b);
    return;
}

// set functions take an already projected x_value
// Private
void Point::set_x(double x_in) {
	x = x_in;
}

// Private
void Point::set_y(double y_in) {
	y = y_in;
}

void Point::set_z(double z_in) {
	z = z_in;
}

void Point::set_R_local(double R_local_in) {
	R_local = R_local_in;
}

double Point::calc_dist_xy(Point in){
    
    return std::sqrt(pow(x - in.get_x(), 2) + pow(y - in.get_y(),2));
    
}


void Point::set_xyz(double xin, double yin, double zin){
    x = xin;
    y = yin;
    z = zin;
}

void Point::set_Point(double gdLatIN, double lonIN, double z_in, double R_local_in) {
    // Save the value if you want
//    gdlat = gdlatIN;
    
    // Calculate the projection
    R_local = R_local_in;
    z = z_in;
//    x = R_equator * lonIN;  // R_equator * longitude_angle
//    y = R_local * gdLatIN;                // Local Radius of Earth * gdlat

    set_xy_from_latlon(gdLatIN, lonIN);
}


void Point::set_xy_from_latlon(double gdLatIN, double lonIN){
    // I think it's not actually necessary to convert the outer degrees into radians because all that matters
    //  is whatever you do you do it consistently.  However, when thinking about the swinging arm length, if
    //  we DO convert these here, then they will give us distances that are easy to interpret i.e. [km]
    
    // Actually, need a projection that preserves area so that each grid cell will have equal area.  This is
    //   required when calculating the exclusion points because we sort the points based on their probability mass.
    //   This sort is only meaningful if all cells have the same area.  FIGURE THIS OUT!!!
    
    // This is an Equidistant Cylindrical Projection
    // It is neither conformal, nor equal-area, and is equidistant only along the reference lines...however, it is fast.
    // http://www.progonos.com/furuti/MapProj/Normal/CartHow/HowER_W12/howER_W12.html
    
    x = refRadius * (lonIN - refLon) * cos(refLat);
    y = refRadius * gdLatIN;
}

//void Point::set_lon(double lonIN) {
//    lon = lonIN;
//}
//
double Point::get_gdLatDeg() {
    // Behind the scenes, everything is radians.
    //  But I know that this is only used at the end, so give it degrees
    return ((y / refRadius) * (180 / PI));
//    return ((180/PI)*(y/R_local));
}

double Point::get_lonDeg() {
//    cout << "refLon = " << refLon *(180/PI) << std::endl;
//    cout << "x/refRadius = " << x *(180/PI) / (refRadius) << std::endl;
//    cout << "refLat and cos() = " << refLat*(180/PI) << "   " << cos(refLat) << std::endl;
//    cout << "value = " << ((refLon + ((x / (refRadius * cos(refLat) )) * (180/PI)) ) -90) << std::endl;
    return ((refLon + ((x / (refRadius * cos(refLat) )) * (180/PI)) ));
//    return ((180/PI) * (x/R_equator));
}

//double Point::get_refLat(){
//    return refLat;
//}
//
//double Point::get_refLon(){
//    return refLon;
//}




// Get functions return an x value
double Point::get_x() {
	return x;
}

double Point::get_y() {
	return y;
}

double Point::get_z() {
	return z;
}

double Point::get_R_local() {
	return R_local;
}


void Point::set_id(int id_in) {
	num_id = id_in;
}

int Point::get_id() {
	return num_id;
}

void Point::set_sortBy(int val){
    // sort by x: val = 1
    // sort by y: val = 2
    // sort by z: val = 3
    // sort by y,x: val = 4
    sortBy = val;
}

//void Point::printMapProperties(){
//    
//    return;
//}
