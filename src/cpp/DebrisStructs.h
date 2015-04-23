/*
 *  DebrisStructs.h
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 7/14/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef DEBRIS_STRUCTS
#define DEBRIS_STRUCTS

struct DebrisCatalogPiece {
	int ID;
	
	int NumLo;
	int NumHi;
	
	double WeightLo;	//Kg
	double WeightHi;
	double WeightTotal;
	
	double RefAreaHi;	//m^2
	double RefAreaLo;
	
	double BallisticLo;		// kg/m^2
	double BallisticHi;
	
	double VelIncrementLo;	// km/s
	double VelIncrementHi;
	
	int CdCode;
	int ClCode;
};

struct DebrisPiece {
	double RefArea;
	int CdCode;
	int ClCode;
	
	//Can't figure out how to get a ColVec object here and I don't want to use a point for fear of memory leaks
	double dVx,dVy,dVz,Mass;
};


#endif