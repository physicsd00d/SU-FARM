/*
 *  Random_Number.cpp
 *  Traj_Prop
 *
 *  Created by Thomas Colvin on 7/1/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "Random_Number.h"
#include <iostream>
using std::endl;
using std::cout;

Random_Number::Random_Number(){
	seed = get_random_seed();
	
//	cout << "~~~~~~~~~~~~~~~USING A FIXED SEED!!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//	seed = 80;
	
	rng = gsl_rng_alloc (gsl_rng_taus);       // allocate the rng
	gsl_rng_set (rng, seed);  // seed the rng
}

Random_Number::~Random_Number(){
	cout << "Freeing the random number generator!!!!!!!!!!!!!!!!!!!!!!" << endl;
	gsl_rng_free (rng);
}

unsigned long int Random_Number::get_random_seed(){
	struct timeval tv;
	unsigned long int seed;       // seed for random number generator
	
	//	cout << "Using time from gettimeofday() for seed . . ." << endl;
	gettimeofday (&tv, 0);
	seed = tv.tv_sec + tv.tv_usec;  // build a pseudo-random seed
	
	//fails at cur_index==65 for seed = 1313439918; arm=50; num_runs= 2000;
	//   NOT related to arm length.  tried arm=200 and still failed.

	
//	seed = 1;
	
	cout << "seed = " << seed << endl;
//	return 1313476854;
	return seed;
}

double Random_Number::generate_number(double sigma) {
	return gsl_ran_gaussian(rng, sigma);
}

double Random_Number::generate_random_gaussian(double mean, double sigma) {
	return mean + gsl_ran_gaussian(rng, sigma);
}

double Random_Number::gaussian_cdf(double x, double sigma){
	return gsl_cdf_gaussian_P(x,sigma);
};

double Random_Number::generate_random_uniform(double lo, double hi) {
	return gsl_ran_flat(rng, lo, hi);
}
