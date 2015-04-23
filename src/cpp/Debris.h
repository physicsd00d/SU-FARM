/*
 *  Debris.h
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 7/13/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DEBRIS_CLASS
#define DEBRIS_CLASS

#define DEBUG_OPTION false

#define INTxx int

#include <vector>
#include <iostream>
using std::cout;
using std::endl;

#include <cmath>

#include "Random_Number.h"
#include "DebrisStructs.h"


//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/matrix.hpp>

//#define ColVec boost::numeric::ublas::vector<double>


class Debris {
	public:
	Debris();
	~Debris();
	
	void LoadShuttleCatalog();

	int GenerateRandomPieces();
	DebrisPiece GetDebrisPiece(int ix);
	int GetNumPieces();

	// Filter out a certain range of debris pieces for debugging purposes
	void ReduceRandomPieces(double ballHi, double ballLo);

	
	int num_types_of_pieces;
	
	std::vector <DebrisCatalogPiece> ShuttleCatalog;
	std::vector <DebrisPiece> RandomPieces;
	
	private:
	Random_Number *rand_num;
	int total_num;
};


#endif DEBRIS_CLASS
