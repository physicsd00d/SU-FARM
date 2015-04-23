/*
 *  Random_Number.h
 *  Traj_Prop
 *
 *  Created by Thomas Colvin on 7/1/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef RANDOM_NUMBER_CLASS
#define RANDOM_NUMBER_CLASS

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include <sys/time.h>           // for the random number seed


class Random_Number{
private:
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	
	unsigned long int seed;
	
	public:
	Random_Number();
	~Random_Number();
	unsigned long int get_random_seed();
	double generate_number(double sigma);
	double generate_random_gaussian(double mean, double sigma);
	double generate_random_uniform(double lo, double hi);

	double gaussian_cdf(double x, double sigma);

};


#endif
