/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "EasyClock.h"
#include "NormalDistribution.h"
#include "BayesianMultivariateRegression.h"
#include "Matrix.h"
#include "Vector.h"
#include "Random.h"
#include <vector>
#include <memory>

// -- Continuous environments -- //
#include "Environment.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "PuddleWorld.h"

// -- Randomness -- //
#include "RandomNumberFile.h"
#include "MersenneTwister.h"

// -- Usual includes -- //
#include <cstring>
#include <getopt.h>

static const char* const help_text = "Usage bayesian multivariate regression [options] \n\
\nOptions: \n\
--environment:	{PuddleWorld, MountainCar, Pendulum}\n\
--randomness:	environment randomness (*0.0) \n\
--n_samples:	Total number of representative samples \n\
--n_episodes:	Total number of episodes \n\
--horizon:		Total number of steps per episode \n\
--b_functions:  Basis function usage or not \n\
--grids:		number of grid intervals \n\
--a:			linear model parameter \n\
--N0:			linear model parameter \n\
\n\
* denotes default parameters \n\
\n";

int main(int argc, char* argv[])
{
	int a			= 0.1;
	int N0			= 0.1;
	
	
	RandomNumberGenerator* rng;
	RandomNumberGenerator* environment_rng;
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;
	
    srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);
    rng->manualSeed(1361690241);
	
	std::cout << "Starting test program" << std::endl;


	int m	= n_states + 1; //Input dimensions (input state dimension plus a dummy state)
	int d_r = 1;			// Reward dimensions
	int d_s = n_states;		//Output states dimensions 	
	
	RBFBasisSet* RBFs = NULL;
	if( b_functions == 1) {
		//Lower and upper environment bounds
		Vector S_L	= environment->StateLowerBound();
		S_L.print(stdout);
		Vector S_U	= environment->StateUpperBound();
		S_U.print(stdout);
		
		std::cout << "Creating Radial basis functions..." << std::endl;
		EvenGrid Discretisation(S_L, S_U, grids);
		RBFs = new RBFBasisSet(Discretisation);

		m = RBFs->size() + 1; // redefinition of the input dimensions (size of the created basis functions plus a dummy state)
		
		std::cout << "# Number of basis functions: " << (RBFs->size() + 1)<< std::endl;
	}
	
	std::cout << "# Creating " << n_actions << " bayesian multivariate models..." << std::endl;
	
	
	//if (argc != 6) {
//		fprintf(stderr, "Usage: Bayesian Multivariate Regression points dimensions a N_0\n");
//		exit(-1);
//	}
//	int n_points = atoi(argv[1]);
//	int n_dimensions_x = atoi(argv[2]);
//	int n_dimensions_y = atoi(argv[3]);
//	real a = atof(argv[4]);
//	real N0 = atof(argv[5]);
//
//
//	BayesianMultivariateRegression bmr(n_dimensions_x, n_dimensions_y, N0 * Matrix::Unity(n_dimensions_y, n_dimensions_y), N0, a);
//	std::auto_ptr<Distribution> distribution_x(new NormalDistribution(0.0,1.0));
//	std::auto_ptr<Distribution> distribution_y(new NormalDistribution(1,0));
//   
//	std::vector<Vector> X(n_points); ///< Input vector
//	std::vector<Vector> Y(n_points); ///< Output vector
//	for( int i = 0; i < n_points; ++i) {
//		X[i].Resize(n_dimensions_x);
//		Y[i].Resize(n_dimensions_y);
//		for(int j = 0; j < n_dimensions_y; ++j) {
//			real ZX  = distribution_x->generate();
//			X[i][j] = urandom(2.0,7.0);
//			Y[i][j] = pow(X[i][j],2.0) + ZX;
//		}
//		X[i][n_dimensions_y] = 1.0;
//	}
//
//	for( int i = 0; i < n_points; ++i) {
//		bmr.AddElement(Y[i],X[i]);
//	}
//	
//	const Matrix W1 = bmr.generate();
//	W1.print(stdout);
//
//	const Matrix W2 = bmr.generate();
//	W2.print(stdout);
//
//	const Matrix W3 = bmr.generate();
//	W3.print(stdout);
//	
//	std::vector<Vector> R1(n_points); 
//	std::vector<Vector> R2(n_points); 
//	std::vector<Vector> R3(n_points); 
//	for(int i = 0; i < n_points; ++i){
//		R1[i] = W1 * X[i];
//		R2[i] = W2 * X[i];
//		R3[i] = W3 * X[i];
//		printf("Input|");
//		X[i].print(stdout);
//		printf("Output|");
//		Y[i].print(stdout);
//		printf("Prediction1|");
//		R1[i].print(stdout);
//		printf("Prediction2|");
//		R2[i].print(stdout);
//		printf("Prediction3|");
//		R3[i].print(stdout);
//	}
}

#endif
