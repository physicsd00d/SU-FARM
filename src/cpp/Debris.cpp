/*
 *  Debris.cpp
 *  Prop3Dof
 *
 *  Created by Thomas Colvin on 7/13/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "Debris.h"

Debris::Debris() {

	//Set up the RNG
	rand_num = new Random_Number();
	total_num = 0;
	
	LoadShuttleCatalog();
	
	return;
}

Debris::~Debris() {
	
	//delete rand_num;
	
	return;
}

void Debris::LoadShuttleCatalog() {
	
	//Empty Piece
	DebrisCatalogPiece EmptyPiece;
	
	num_types_of_pieces = 25;
	ShuttleCatalog.assign(num_types_of_pieces, EmptyPiece);
	

	num_types_of_pieces = 25;
	ShuttleCatalog[0].NumLo = 1;
	ShuttleCatalog[0].NumHi = 1;
	ShuttleCatalog[0].WeightLo = 2948.348;
	ShuttleCatalog[0].WeightHi = 2948.348;
	ShuttleCatalog[0].WeightTotal = 2948.348;
	ShuttleCatalog[0].RefAreaLo = 55.18440576;
	ShuttleCatalog[0].RefAreaHi = 55.18440576;
	ShuttleCatalog[0].BallisticLo = 39.0593892299;
	ShuttleCatalog[0].BallisticHi = 78.1187784598;
	ShuttleCatalog[0].VelIncrementLo = 0.001524;
	ShuttleCatalog[0].VelIncrementHi = 0.001524;
	ShuttleCatalog[0].CdCode = 1;
	ShuttleCatalog[0].ClCode = 1;
	
	ShuttleCatalog[1].NumLo = 1;
	ShuttleCatalog[1].NumHi = 1;
	ShuttleCatalog[1].WeightLo = 19232.3008;
	ShuttleCatalog[1].WeightHi = 19232.3008;
	ShuttleCatalog[1].WeightTotal = 192323.008;
	ShuttleCatalog[1].RefAreaLo = 55.18440576;
	ShuttleCatalog[1].RefAreaHi = 55.18440576;
	ShuttleCatalog[1].BallisticLo = 73.236354806;
	ShuttleCatalog[1].BallisticHi = 209.944217111;
	ShuttleCatalog[1].VelIncrementLo = 0.0033528;
	ShuttleCatalog[1].VelIncrementHi = 0.0033528;
	ShuttleCatalog[1].CdCode = 2;
	ShuttleCatalog[1].ClCode = 2;
	
	ShuttleCatalog[2].NumLo = 1;
	ShuttleCatalog[2].NumHi = 1;
	ShuttleCatalog[2].WeightLo = 9253.2768;
	ShuttleCatalog[2].WeightHi = 9253.2768;
	ShuttleCatalog[2].WeightTotal = 9253.2768;
	ShuttleCatalog[2].RefAreaLo = 55.18440576;
	ShuttleCatalog[2].RefAreaHi = 55.18440576;
	ShuttleCatalog[2].BallisticLo = 68.3539311523;
	ShuttleCatalog[2].BallisticHi = 292.945419224;
	ShuttleCatalog[2].VelIncrementLo = 0.0112776;
	ShuttleCatalog[2].VelIncrementHi = 0.0112776;
	ShuttleCatalog[2].CdCode = 3;
	ShuttleCatalog[2].ClCode = 3;
	
	ShuttleCatalog[3].NumLo = 1;
	ShuttleCatalog[3].NumHi = 1;
	ShuttleCatalog[3].WeightLo = 240.40376;
	ShuttleCatalog[3].WeightHi = 240.40376;
	ShuttleCatalog[3].WeightTotal = 240.40376;
	ShuttleCatalog[3].RefAreaLo = 1.11483648;
	ShuttleCatalog[3].RefAreaHi = 1.11483648;
	ShuttleCatalog[3].BallisticLo = 268.533300955;
	ShuttleCatalog[3].BallisticHi = 927.66049421;
	ShuttleCatalog[3].VelIncrementLo = 0.064008;
	ShuttleCatalog[3].VelIncrementHi = 0.064008;
	ShuttleCatalog[3].CdCode = 4;
	ShuttleCatalog[3].ClCode = 0;
	
	ShuttleCatalog[4].NumLo = 2;
	ShuttleCatalog[4].NumHi = 2;
	ShuttleCatalog[4].WeightLo = 61.23492;
	ShuttleCatalog[4].WeightHi = 147.4174;
	ShuttleCatalog[4].WeightTotal = 208.65232;
	ShuttleCatalog[4].RefAreaLo = 1.11483648;
	ShuttleCatalog[4].RefAreaHi = 2.7870912;
	ShuttleCatalog[4].BallisticLo = 68.3539311523;
	ShuttleCatalog[4].BallisticHi = 263.650877302;
	ShuttleCatalog[4].VelIncrementLo = 0.064008;
	ShuttleCatalog[4].VelIncrementHi = 0.064008;
	ShuttleCatalog[4].CdCode = 4;
	ShuttleCatalog[4].ClCode = 0;
	
	ShuttleCatalog[5].NumLo = 3;
	ShuttleCatalog[5].NumHi = 3;
	ShuttleCatalog[5].WeightLo = 9.07184;
	ShuttleCatalog[5].WeightHi = 18.14368;
	ShuttleCatalog[5].WeightTotal = 45.3592;
	ShuttleCatalog[5].RefAreaLo = 0.27870912;
	ShuttleCatalog[5].RefAreaHi = 0.65032128;
	ShuttleCatalog[5].BallisticLo = 34.1769655762;
	ShuttleCatalog[5].BallisticHi = 112.295744036;
	ShuttleCatalog[5].VelIncrementLo = 0.064008;
	ShuttleCatalog[5].VelIncrementHi = 0.064008;
	ShuttleCatalog[5].CdCode = 4;
	ShuttleCatalog[5].ClCode = 0;
	
	ShuttleCatalog[6].NumLo = 3;
	ShuttleCatalog[6].NumHi = 3;
	ShuttleCatalog[6].WeightLo = 3.628736;
	ShuttleCatalog[6].WeightHi = 6.80388;
	ShuttleCatalog[6].WeightTotal = 14.968536;
	ShuttleCatalog[6].RefAreaLo = 0.18580608;
	ShuttleCatalog[6].RefAreaHi = 0.32516064;
	ShuttleCatalog[6].BallisticLo = 24.4121182687;
	ShuttleCatalog[6].BallisticHi = 78.1187784598;
	ShuttleCatalog[6].VelIncrementLo = 0.01524;
	ShuttleCatalog[6].VelIncrementHi = 0.064008;
	ShuttleCatalog[6].CdCode = 4;
	ShuttleCatalog[6].ClCode = 0;
	
	ShuttleCatalog[7].NumLo = 2;
	ShuttleCatalog[7].NumHi = 2;
	ShuttleCatalog[7].WeightLo = 4.53592;
	ShuttleCatalog[7].WeightHi = 6.80388;
	ShuttleCatalog[7].WeightTotal = 11.3398;
	ShuttleCatalog[7].RefAreaLo = 0.2322576;
	ShuttleCatalog[7].RefAreaHi = 0.32516064;
	ShuttleCatalog[7].BallisticLo = 24.4121182687;
	ShuttleCatalog[7].BallisticHi = 78.1187784598;
	ShuttleCatalog[7].VelIncrementLo = 0.064008;
	ShuttleCatalog[7].VelIncrementHi = 0.064008;
	ShuttleCatalog[7].CdCode = 4;
	ShuttleCatalog[7].ClCode = 0;
	
	ShuttleCatalog[8].NumLo = 2;
	ShuttleCatalog[8].NumHi = 2;
	ShuttleCatalog[8].WeightLo = 0.3628736;
	ShuttleCatalog[8].WeightHi = 0.3628736;
	ShuttleCatalog[8].WeightTotal = 0.7257472;
	ShuttleCatalog[8].RefAreaLo = 0.018580608;
	ShuttleCatalog[8].RefAreaHi = 0.018580608;
	ShuttleCatalog[8].BallisticLo = 19.5296946149;
	ShuttleCatalog[8].BallisticHi = 146.472709612;
	ShuttleCatalog[8].VelIncrementLo = 0.01524;
	ShuttleCatalog[8].VelIncrementHi = 0.01524;
	ShuttleCatalog[8].CdCode = 5;
	ShuttleCatalog[8].ClCode = 0;
	
	ShuttleCatalog[9].NumLo = 2;
	ShuttleCatalog[9].NumHi = 2;
	ShuttleCatalog[9].WeightLo = 0.3628736;
	ShuttleCatalog[9].WeightHi = 0.3628736;
	ShuttleCatalog[9].WeightTotal = 0.7257472;
	ShuttleCatalog[9].RefAreaLo = 0.018580608;
	ShuttleCatalog[9].RefAreaHi = 0.018580608;
	ShuttleCatalog[9].BallisticLo = 19.5296946149;
	ShuttleCatalog[9].BallisticHi = 146.472709612;
	ShuttleCatalog[9].VelIncrementLo = 0.01524;
	ShuttleCatalog[9].VelIncrementHi = 0.01524;
	ShuttleCatalog[9].CdCode = 5;
	ShuttleCatalog[9].ClCode = 0;
	
	ShuttleCatalog[10].NumLo = 4;
	ShuttleCatalog[10].NumHi = 4;
	ShuttleCatalog[10].WeightLo = 0.5896696;
	ShuttleCatalog[10].WeightHi = 0.5896696;
	ShuttleCatalog[10].WeightTotal = 2.3586784;
	ShuttleCatalog[10].RefAreaLo = 0.027870912;
	ShuttleCatalog[10].RefAreaHi = 0.027870912;
	ShuttleCatalog[10].BallisticLo = 19.5296946149;
	ShuttleCatalog[10].BallisticHi = 146.472709612;
	ShuttleCatalog[10].VelIncrementLo = 0.064008;
	ShuttleCatalog[10].VelIncrementHi = 0.064008;
	ShuttleCatalog[10].CdCode = 5;
	ShuttleCatalog[10].ClCode = 0;
	
	ShuttleCatalog[11].NumLo = 1;
	ShuttleCatalog[11].NumHi = 1;
	ShuttleCatalog[11].WeightLo = 20.41164;
	ShuttleCatalog[11].WeightHi = 20.41164;
	ShuttleCatalog[11].WeightTotal = 20.41164;
	ShuttleCatalog[11].RefAreaLo = 1.3935456;
	ShuttleCatalog[11].RefAreaHi = 1.3935456;
	ShuttleCatalog[11].BallisticLo = 19.5296946149;
	ShuttleCatalog[11].BallisticHi = 19.5296946149;
	ShuttleCatalog[11].VelIncrementLo = 0.01524;
	ShuttleCatalog[11].VelIncrementHi = 0.01524;
	ShuttleCatalog[11].CdCode = 6;
	ShuttleCatalog[11].ClCode = 0;
	
	ShuttleCatalog[12].NumLo = 1;
	ShuttleCatalog[12].NumHi = 1;
	ShuttleCatalog[12].WeightLo = 188.24068;
	ShuttleCatalog[12].WeightHi = 188.24068;
	ShuttleCatalog[12].WeightTotal = 188.24068;
	ShuttleCatalog[12].RefAreaLo = 13.0064256;
	ShuttleCatalog[12].RefAreaHi = 13.0064256;
	ShuttleCatalog[12].BallisticLo = 19.5296946149;
	ShuttleCatalog[12].BallisticHi = 19.5296946149;
	ShuttleCatalog[12].VelIncrementLo = 0.01524;
	ShuttleCatalog[12].VelIncrementHi = 0.01524;
	ShuttleCatalog[12].CdCode = 6;
	ShuttleCatalog[12].ClCode = 0;
	
	ShuttleCatalog[13].NumLo = 1;
	ShuttleCatalog[13].NumHi = 1;
	ShuttleCatalog[13].WeightLo = 40.82328;
	ShuttleCatalog[13].WeightHi = 40.82328;
	ShuttleCatalog[13].WeightTotal = 40.82328;
	ShuttleCatalog[13].RefAreaLo = 2.7870912;
	ShuttleCatalog[13].RefAreaHi = 2.7870912;
	ShuttleCatalog[13].BallisticLo = 19.5296946149;
	ShuttleCatalog[13].BallisticHi = 19.5296946149;
	ShuttleCatalog[13].VelIncrementLo = 0.064008;
	ShuttleCatalog[13].VelIncrementHi = 0.064008;
	ShuttleCatalog[13].CdCode = 6;
	ShuttleCatalog[13].ClCode = 0;
	
	ShuttleCatalog[14].NumLo = 1;
	ShuttleCatalog[14].NumHi = 1;
	ShuttleCatalog[14].WeightLo = 331.12216;
	ShuttleCatalog[14].WeightHi = 331.12216;
	ShuttleCatalog[14].WeightTotal = 331.12216;
	ShuttleCatalog[14].RefAreaLo = 25.0838208;
	ShuttleCatalog[14].RefAreaHi = 25.0838208;
	ShuttleCatalog[14].BallisticLo = 19.5296946149;
	ShuttleCatalog[14].BallisticHi = 19.5296946149;
	ShuttleCatalog[14].VelIncrementLo = 0.064008;
	ShuttleCatalog[14].VelIncrementHi = 0.064008;
	ShuttleCatalog[14].CdCode = 6;
	ShuttleCatalog[14].ClCode = 0;
	
	ShuttleCatalog[15].NumLo = 1;
	ShuttleCatalog[15].NumHi = 1;
	ShuttleCatalog[15].WeightLo = 390.08912;
	ShuttleCatalog[15].WeightHi = 390.08912;
	ShuttleCatalog[15].WeightTotal = 390.08912;
	ShuttleCatalog[15].RefAreaLo = 29.7289728;
	ShuttleCatalog[15].RefAreaHi = 29.7289728;
	ShuttleCatalog[15].BallisticLo = 19.5296946149;
	ShuttleCatalog[15].BallisticHi = 19.5296946149;
	ShuttleCatalog[15].VelIncrementLo = 0.064008;
	ShuttleCatalog[15].VelIncrementHi = 0.064008;
	ShuttleCatalog[15].CdCode = 6;
	ShuttleCatalog[15].ClCode = 0;
	
	ShuttleCatalog[16].NumLo = 8;
	ShuttleCatalog[16].NumHi = 8;
	ShuttleCatalog[16].WeightLo = 0.0453592;
	ShuttleCatalog[16].WeightHi = 2.26796;
	ShuttleCatalog[16].WeightTotal = 9.07184;
	ShuttleCatalog[16].RefAreaLo = 0.0037161216;
	ShuttleCatalog[16].RefAreaHi = 0.09290304;
	ShuttleCatalog[16].BallisticLo = 9.76484730747;
	ShuttleCatalog[16].BallisticHi = 156.23755692;
	ShuttleCatalog[16].VelIncrementLo = 0.064008;
	ShuttleCatalog[16].VelIncrementHi = 0.064008;
	ShuttleCatalog[16].CdCode = 7;
	ShuttleCatalog[16].ClCode = 0;
	
	ShuttleCatalog[17].NumLo = 14;
	ShuttleCatalog[17].NumHi = 14;
	ShuttleCatalog[17].WeightLo = 1.814368;
	ShuttleCatalog[17].WeightHi = 9.979024;
	ShuttleCatalog[17].WeightTotal = 44.452016;
	ShuttleCatalog[17].RefAreaLo = 0.18580608;
	ShuttleCatalog[17].RefAreaHi = 0.83612736;
	ShuttleCatalog[17].BallisticLo = 7.81187784598;
	ShuttleCatalog[17].BallisticHi = 14.6472709612;
	ShuttleCatalog[17].VelIncrementLo = 0.064008;
	ShuttleCatalog[17].VelIncrementHi = 0.064008;
	ShuttleCatalog[17].CdCode = 8;
	ShuttleCatalog[17].ClCode = 0;
	
	ShuttleCatalog[18].NumLo = 16;
	ShuttleCatalog[18].NumHi = 48;
	ShuttleCatalog[18].WeightLo = 0.453592;
	ShuttleCatalog[18].WeightHi = 1.360776;
	ShuttleCatalog[18].WeightTotal = 21.772416;
	ShuttleCatalog[18].RefAreaLo = 0.018580608;
	ShuttleCatalog[18].RefAreaHi = 0.065032128;
	ShuttleCatalog[18].BallisticLo = 16.6002404227;
	ShuttleCatalog[18].BallisticHi = 29.2945419224;
	ShuttleCatalog[18].VelIncrementLo = 0.01524;
	ShuttleCatalog[18].VelIncrementHi = 0.01524;
	ShuttleCatalog[18].CdCode = 8;
	ShuttleCatalog[18].ClCode = 0;
	
	ShuttleCatalog[19].NumLo = 4;
	ShuttleCatalog[19].NumHi = 48;
	ShuttleCatalog[19].WeightLo = 1.587572;
	ShuttleCatalog[19].WeightHi = 19.050864;
	ShuttleCatalog[19].WeightTotal = 77.11064;
	ShuttleCatalog[19].RefAreaLo = 0.4645152;
	ShuttleCatalog[19].RefAreaHi = 5.5741824;
	ShuttleCatalog[19].BallisticLo = 4.39418128836;
	ShuttleCatalog[19].BallisticHi = 5.85890838448;
	ShuttleCatalog[19].VelIncrementLo = 0.0027432;
	ShuttleCatalog[19].VelIncrementHi = 0.0027432;
	ShuttleCatalog[19].CdCode = 9;
	ShuttleCatalog[19].ClCode = 0;
	
	ShuttleCatalog[20].NumLo = 160;
	ShuttleCatalog[20].NumHi = 160;
	ShuttleCatalog[20].WeightLo = 0.0453592;
	ShuttleCatalog[20].WeightHi = 0.0907184;
	ShuttleCatalog[20].WeightTotal = 10.886208;
	ShuttleCatalog[20].RefAreaLo = 0.018580608;
	ShuttleCatalog[20].RefAreaHi = 0.027870912;
	ShuttleCatalog[20].BallisticLo = 1.95296946149;
	ShuttleCatalog[20].BallisticHi = 3.41769655762;
	ShuttleCatalog[20].VelIncrementLo = 0.01524;
	ShuttleCatalog[20].VelIncrementHi = 0.01524;
	ShuttleCatalog[20].CdCode = 8;
	ShuttleCatalog[20].ClCode = 0;
	
	ShuttleCatalog[21].NumLo = 16;
	ShuttleCatalog[21].NumHi = 32;
	ShuttleCatalog[21].WeightLo = 0.680388;
	ShuttleCatalog[21].WeightHi = 1.13398;
	ShuttleCatalog[21].WeightTotal = 25.401152;
	ShuttleCatalog[21].RefAreaLo = 0.501676416;
	ShuttleCatalog[21].RefAreaHi = 0.501676416;
	ShuttleCatalog[21].BallisticLo = 1.46472709612;
	ShuttleCatalog[21].BallisticHi = 3.90593892299;
	ShuttleCatalog[21].VelIncrementLo = 0.01524;
	ShuttleCatalog[21].VelIncrementHi = 0.01524;
	ShuttleCatalog[21].CdCode = 9;
	ShuttleCatalog[21].ClCode = 0;
	
	ShuttleCatalog[22].NumLo = 50;
	ShuttleCatalog[22].NumHi = 50;
	ShuttleCatalog[22].WeightLo = 0.4082328;
	ShuttleCatalog[22].WeightHi = 8.618248;
	ShuttleCatalog[22].WeightTotal = 75.749864;
	ShuttleCatalog[22].RefAreaLo = 0.09290304;
	ShuttleCatalog[22].RefAreaHi = 0.65032128;
	ShuttleCatalog[22].BallisticLo = 244.121182687;
	ShuttleCatalog[22].BallisticHi = 976.484730747;
	ShuttleCatalog[22].VelIncrementLo = 0.01524;
	ShuttleCatalog[22].VelIncrementHi = 0.064008;
	ShuttleCatalog[22].CdCode = 8;
	ShuttleCatalog[22].ClCode = 0;
	
	ShuttleCatalog[23].NumLo = 50;
	ShuttleCatalog[23].NumHi = 200;
	ShuttleCatalog[23].WeightLo = 0.002721552;
	ShuttleCatalog[23].WeightHi = 0.0226796;
	ShuttleCatalog[23].WeightTotal = 0.9979024;
	ShuttleCatalog[23].RefAreaLo = 0.0009290304;
	ShuttleCatalog[23].RefAreaHi = 0.0074322432;
	ShuttleCatalog[23].BallisticLo = 2.44121182687;
	ShuttleCatalog[23].BallisticHi = 3.90593892299;
	ShuttleCatalog[23].VelIncrementLo = 0.70104;
	ShuttleCatalog[23].VelIncrementHi = 0.70104;
	ShuttleCatalog[23].CdCode = 8;
	ShuttleCatalog[23].ClCode = 0;
	
	ShuttleCatalog[24].NumLo = 100;
	ShuttleCatalog[24].NumHi = 500;
	ShuttleCatalog[24].WeightLo = 0.00453592;
	ShuttleCatalog[24].WeightHi = 2.26796;
	ShuttleCatalog[24].WeightTotal = 453.592;
	ShuttleCatalog[24].RefAreaLo = 9.290304e-05;
	ShuttleCatalog[24].RefAreaHi = 0.027870912;
	ShuttleCatalog[24].BallisticLo = 4.88242365374;
	ShuttleCatalog[24].BallisticHi = 97.6484730747;
	ShuttleCatalog[24].VelIncrementLo = 0.01524;
	ShuttleCatalog[24].VelIncrementHi = 0.064008;
	ShuttleCatalog[24].CdCode = 1;
	ShuttleCatalog[24].ClCode = 0;

	return;
}

int Debris::GenerateRandomPieces(){

	int *num_of_each_type = new int [num_types_of_pieces];
	total_num = 0;
	
	// First determine how many pieces of debris we will have
	for (int i = 0; i < num_types_of_pieces; i++) {
		double Lo = ShuttleCatalog[i].NumLo;
		double Hi = ShuttleCatalog[i].NumHi;
		
		if (DEBUG_OPTION) {
			num_of_each_type[i] = 5; }
		else {
			num_of_each_type[i] = (INTxx) rand_num->generate_random_uniform(Lo,Hi); }

		total_num += num_of_each_type[i]; }
	
	cout << "%Total num = " << total_num << endl;
	
	// Allocate a vector of that length
	DebrisPiece EmptyPiece;
	RandomPieces.assign(total_num,EmptyPiece);
	
	// Populate the vector
	int cur_piece = 0;
	for (int ix = 0; ix < num_types_of_pieces; ix++) {
		
		for (int j = 0; j < num_of_each_type[ix]; j++) {
			// Not worrying about hitting the 'total mass' requirement
			RandomPieces[cur_piece].Mass  
				= rand_num->generate_random_uniform(ShuttleCatalog[ix].WeightLo, ShuttleCatalog[ix].WeightHi);
			
			// Reference Area
			RandomPieces[cur_piece].RefArea
				= rand_num->generate_random_uniform(ShuttleCatalog[ix].RefAreaLo, ShuttleCatalog[ix].RefAreaHi);
			
			// Record which Lift/Drag code to use
			RandomPieces[cur_piece].CdCode = ShuttleCatalog[ix].CdCode;
			RandomPieces[cur_piece].ClCode = ShuttleCatalog[ix].ClCode;
			
			// Generate a random unit vector
			double ux = rand_num->generate_random_uniform(0, 1);
			double uy = rand_num->generate_random_uniform(0, 1);
			double uz = rand_num->generate_random_uniform(0, 1);
			double norm = std::sqrt(ux*ux + uy*uy + uz*uz);
			
			// Compute the velocity increment
			RandomPieces[cur_piece].dVx
				= rand_num->generate_random_uniform(ShuttleCatalog[ix].VelIncrementLo, ShuttleCatalog[ix].VelIncrementHi) * ux/norm;
			RandomPieces[cur_piece].dVy
				= rand_num->generate_random_uniform(ShuttleCatalog[ix].VelIncrementLo, ShuttleCatalog[ix].VelIncrementHi) * uy/norm;
			RandomPieces[cur_piece].dVz
				= rand_num->generate_random_uniform(ShuttleCatalog[ix].VelIncrementLo, ShuttleCatalog[ix].VelIncrementHi) * uz/norm;
			
			cur_piece++;
		} }
	
	return total_num;
}

DebrisPiece Debris::GetDebrisPiece(int ix) {
	
	if ( (ix < total_num) && (ix >=0) ) {
		return RandomPieces[ix]; }
	else {
		std::cerr << "This piece doesn't exist!!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		exit(2352345);	}
	
}

int Debris::GetNumPieces(){
	
	return total_num;
}

void Debris::ReduceRandomPieces(double ballHi, double ballLo){
	// Ballistic Coefficient = M/(Cd * ARef)
	// Assuming Cd = 0.5 for now

	double Cd = 0.5;
	std::vector <DebrisPiece> TempPieces;

	DebrisPiece CurrentPiece;
	for (int ix = 0; ix < total_num; ix++){
		CurrentPiece = GetDebrisPiece(ix);

		double ballcoeff = CurrentPiece.Mass/(CurrentPiece.RefArea * Cd);

		if ((ballcoeff <= ballHi) && (ballcoeff >= ballLo)){
			TempPieces.push_back(CurrentPiece);
		}

	}

	// Replace the old vector with the abridged pieces
	RandomPieces.clear();
	RandomPieces = TempPieces;
	total_num = (INTxx) RandomPieces.size();

	return;
}






















