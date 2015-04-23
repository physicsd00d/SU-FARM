/*
 *  Point.h
 *  Traj_Prop
 *
 *  Created by Thomas Colvin on 8/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef POINT_CLASS
#define POINT_CLASS

#define R_equator 6378.1370   //[km]
#define R_polar 6356.7523	  //[km]
#define PI 3.141592653589793


#include <iostream>
using std::ostream;
using std::cout;

#include <cmath>


class Point{
	
private:
	double x, y, z, R_local;		//[km] Adding R_local changes the size of the object, which MAY mess up the reading in functions
//    double gdlat, lon;              //[rad] The lat and lon of a point
    
    //[rad] The center point of the projection
    const static double refLat = 38 * (PI/180);     // This is basically the Lat at which we're going to approximate the radius of the earth
    const static double refLon = 0.;    // You really want this to be zero because we're measuring our Longitude from Greenwich
    
    // The radius at the center point
    double refRadius;
    
	int num_id;

    // User is not currently allowed to directly manipulate the projected data
	void set_x(double x_in);
	void set_y(double y_in);
    void calcRefRadius();
	
public:
	Point(double gdLatIN, double lonIN, double z_in, double R_local_in);    //(rad, rad, km, km)
	Point(const Point &object);
    Point();
	//	~Point();
	
    void set_Point(double gdLatIN, double lonIN, double z_in, double R_local_in);


	void set_z(double z_in);
	void set_R_local(double R_local_in);

	double get_x();
	double get_y();
	double get_z();
	double get_R_local();
    
    double calc_dist_xy(Point obj);
    
    void set_gdlat(double gdlatIN);
    void set_lon(double lonIN);
    
//    void set_x_and_y(double gdLatIN, double lonIN);
    void set_xy_from_latlon(double gdLatIN, double lonIN);
    void set_xyz(double xin, double yin, double zin);

    
    double get_gdLatDeg();
    double get_lonDeg();
    double get_refLat();
    double get_refLon();
    
    double calcAzimuth(const Point &ptFinal);
    void rotateThisAboutAnotherPoint(const Point &refPoint, double rotAngleRad);


	void set_id(int id_in);
	int get_id();
    
    static int sortBy;
    void set_sortBy(int val);
    
    // Extra information about the debris at this point
    double Vx_km_s;
    double Vy_km_s;
    double Vz_km_s;
    double mass_kg;
    double area_km2;

	
	// Operators
	friend ostream& operator<<(ostream& output, const Point& p);
	Point operator+ (Point rhs);
	Point operator- (Point rhs);
	bool operator== (Point rhs);
	bool operator < (Point rhs);

	
	
	//	bool myfunction (Point p1,Point p2);
	
};









#endif