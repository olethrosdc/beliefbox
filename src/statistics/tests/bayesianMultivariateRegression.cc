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
--data:         Data file\n\
--horizon:		Total number of steps per episode \n\
--b_functions:  Basis function usage or not \n\
--grids:		number of grid intervals \n\
--a:			linear model parameter (0.1)\n\
--N0:			linear model parameter (0.1)\n\
\n\
* denotes default parameters \n\
\n";

int main(int argc, char* argv[])
{
	
	RandomNumberGenerator* rng;
    MersenneTwisterRNG mersenne_twister;
    rng = (RandomNumberGenerator*) &mersenne_twister;
	
    srand48(34987235);
    srand(34987235);
    setRandomSeed(34987235);
    rng->manualSeed(1361690241);
	
	std::cout << "# Starting test program" << std::endl;


	if (argc != 6) {
		fprintf(stderr, "Usage: bayesianMultivariateRegression points x_dim y_dim a N_0\n");
		exit(-1);
	}
	int n_points = atoi(argv[1]);
	int n_dimensions_x = atoi(argv[2]);
	int n_dimensions_y = atoi(argv[3]);
	real a = atof(argv[4]);
	real N0 = atof(argv[5]);


	int n_train = (int) ceil(n_points / 2);

	BayesianMultivariateRegression bmr(n_dimensions_x, n_dimensions_y, N0 * Matrix::Unity(n_dimensions_y, n_dimensions_y), N0, a);
	NormalDistribution generator(0.0,1.0);
	Matrix A(n_dimensions_x, n_dimensions_y); // mean 
	Matrix P(n_dimensions_y, n_dimensions_y); // precision
	for(int i = 0; i < n_dimensions_x; ++i) {
		for(int j = 0; j < n_dimensions_y; ++j) {
			A(i,j) = generator.generate();
			P(i,j) = generator.generate();
		}
	}
	P = P * Transpose(P);
	MultivariateNormal noise(Vector::Null(n_dimensions_x), P);
	std::vector<Vector> X(n_points); ///< Input vector
	std::vector<Vector> Y(n_points); ///< Output vector
	for( int i = 0; i < n_points; ++i) {
		X[i].Resize(n_dimensions_x);
		for(int j = 0; j < n_dimensions_x; ++j) {
			X[i](j) = urandom(-1.0,5.0);
		}
		Y[i].Resize(n_dimensions_y);
		Y[i] = A * X[i] + noise.generate();
		printf("# X: ");
		X[i].print(stdout);
		printf("# Y: ");
		Y[i].print(stdout);
	}

	for( int i = 0; i < n_train; ++i) {
		bmr.AddElement(Y[i],X[i]);
	}
	
	const Matrix W1 = bmr.generate();
	W1.print(stdout);
	
	const Matrix W2 = bmr.generate();
	W2.print(stdout);
	
	const Matrix W3 = bmr.generate();
	W3.print(stdout);
	
	std::vector<Vector> R1(n_points); 
	std::vector<Vector> R2(n_points); 
	std::vector<Vector> R3(n_points);
	printf("# Predictions on test points\n");
	printf("#===========================\n");
	for(int i = n_train+1; i < n_points; ++i){
		R1[i] = W1 * X[i];
		R2[i] = W2 * X[i];
		R3[i] = W3 * X[i];
		printf("# Prediction1:");
		R1[i].print(stdout);
		if (0) {
			printf("Input|");
			X[i].print(stdout);
			printf("Output|");
			Y[i].print(stdout);
			printf("# Prediction1:");
			R1[i].print(stdout);
			printf("Prediction2|");
			R2[i].print(stdout);
			printf("Prediction3|");
			R3[i].print(stdout);
		}
		printf("# err: %f\n", EuclideanNorm(R1[i], Y[i]));
	}

	printf("A:\n");
	A.print(stdout);
	printf("W1:\n");
	W1.print(stdout);
	printf("W2:\n");
	W2.print(stdout);
	printf("W3:\n");
	W3.print(stdout);
}

#endif
