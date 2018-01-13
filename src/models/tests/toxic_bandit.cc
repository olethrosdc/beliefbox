// copyright (c) 2015 by Hannes Eriksson <hannese18@gmail.com>
/***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************/

#include "Vector.h"
#include "GaussianProcessBandit.h"
#include "GaussianProcessRPBandit.h"
#include "GaussianProcessTSBandit.h"
#include "LinearBanditTS.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "CompressSensing.h"
#include "MultivariateNormal.h"

// test of GP-UCB and LinearBandit-TS

int readIntoTupleVector(std::vector<std::vector<int> > &descriptorCount, std::vector<std::vector<int> > &descriptorIndex, 
	std::vector<real> &rewardData, const char* fname) {
		std::ifstream input(fname);
		int i = 0;
		for(std::string line; getline(input,line);) {
			int index = line.find(" ");
			real reward = ::atof(line.substr(0,index).c_str());
			std::string tmp = line.substr(index+1,line.length());
			int j = 0;
			while(tmp.length() > 4) {
				int nextIndex = tmp.find(" ");
				if(nextIndex == -1)
					nextIndex = tmp.length()-1;
				std::string nextPair = tmp.substr(0,nextIndex);
				int delimIndex = nextPair.find(":");
				int descIndex = ::atoi(nextPair.substr(0,delimIndex).c_str());
				int descCount = ::atoi(nextPair.substr(delimIndex+1,nextPair.length()).c_str());
				descriptorCount[i][j] = descCount;
				descriptorIndex[i][j] = descIndex;
				tmp = tmp.substr(nextIndex+1,tmp.length());
				j++;
			}
			descriptorCount[i].resize(j);
			descriptorIndex[i].resize(j);
			rewardData[i] = reward;
			i++;
		}
		rewardData.resize(i);
		descriptorCount.resize(i);
		descriptorIndex.resize(i);
		return i;
}

int main(int argc, char* argv[]) {
	std::cout << "Starting\n";
	const int NUM_MOLECULES = 2000;
	const int NUM_DESCRIPTORS = 6000;
	std::vector<std::vector<int> > descriptorCount(NUM_MOLECULES);
	std::vector<std::vector<int> > descriptorIndex(NUM_MOLECULES);
	for(int i=0;i<NUM_MOLECULES;i++) {
		descriptorCount[i].resize(NUM_DESCRIPTORS);
		descriptorIndex[i].resize(NUM_DESCRIPTORS);
	}
	std::vector<real> rewardData(NUM_MOLECULES);
	readIntoTupleVector(descriptorCount,descriptorIndex,rewardData,argv[1]);	

	compressor c = createCompressionMatrix(descriptorIndex);
	std::vector<std::vector<real> > compressed(rewardData.size());

	// COMPRESSED SENSING
	for(int i=0;i<compressed.size();i++) {
		compressed[i] = compressSparseArray(descriptorCount[i],descriptorIndex[i],c);
	}
	Vector Y = Vector(rewardData);
	Matrix X = Matrix(compressed.size(),compressed[0].size());
	for(int i=0;i<compressed.size();i++) {
		X.setRow(i,Vector(compressed[i]));
	}

	/*
	Matrix X = Matrix(25,3);
	Vector Y = Vector(25);
	X(0,0) = 0; X(0,1) = 0; X(0,2) = 0; Y(0) = 8;
	X(1,0) = 1; X(1,1) = 0; X(1,2) = 0; Y(1) = 7;	 
	X(2,0) = 0; X(2,1) = 1; X(2,2) = 0; Y(2) = 7;	
	X(3,0) = 0; X(3,1) = 0; X(3,2) = 1; Y(3) = 7;	 
	X(4,0) = 1; X(4,1) = 1; X(4,2) = 0; Y(4) = 6.5;	 
	X(5,0) = 1; X(5,1) = 0; X(5,2) = 1; Y(5) = 6.5;
	X(6,0) = 0; X(6,1) = 1; X(6,2) = 1; Y(6) = 6.5;
	X(7,0) = 1; X(7,1) = 1; X(7,2) = 1; Y(7) = 6;	 
	X(8,0) = 0.25; X(8,1) = 0.25; X(8,2) = 0; Y(8) = 7.5;	
	X(9,0) = 0; X(9,1) = 0.25; X(9,2) = 0.25; Y(9) = 7.5;	 
	X(10,0) = 0.25; X(10,1) = 0; X(10,2) = 0.25; Y(10) = 7.5;	 
	X(11,0) = 0.5; X(11,1) = 0.5; X(11,2) = 0; Y(11) = 7.1;
	X(12,0) = 0.5; X(12,1) = 0; X(12,2) = 0.5; Y(12) = 7.1;
	X(13,0) = 0; X(13,1) = 0.5; X(13,2) = 0.5; Y(13) = 7.1;	 
	X(14,0) = 0.75; X(14,1) = 0; X(14,2) = 0; Y(14) = 7.25;	
	X(15,0) = 0; X(15,1) = 0.75; X(15,2) = 0; Y(15) = 7.25;	 
	X(16,0) = 0; X(16,1) = 0; X(16,2) = 0.75; Y(16) = 7.25;	 
	X(17,0) = 2; X(17,1) = 0; X(17,2) = 0; Y(17) = 5;
	X(18,0) = 0; X(18,1) = 2; X(18,2) = 0; Y(18) = 5;
	X(19,0) = 0; X(19,1) = 0; X(19,2) = 2; Y(19) = 5;
	X(20,0) = 2; X(20,1) = 2; X(20,2) = 0; Y(20) = 3;	 
	X(21,0) = 2; X(21,1) = 0; X(21,2) = 2; Y(21) = 3;	
	X(22,0) = 0; X(22,1) = 2; X(22,2) = 2; Y(22) = 3;	 
	X(23,0) = 2; X(23,1) = 2; X(23,2) = 2; Y(23) = 2.5;	 
	X(24,0) = 0; X(24,1) = 3; X(24,2) = 0; Y(24) = 4;
	*/
	real scale_length = 2 * std::pow(10.0,2);
	Vector gp_scale_length = Vector(X.Columns()) + scale_length;
	Vector gp_scale_length_dic = Vector(X.Columns()) + scale_length;
	real gp_noise_var = 1.0;
	real gp_sig_var = 1.0;
	real gp_threshold = 0.01;
	real delta = 0.5;
	int iterations = floor(Y.Size() * 0.25);
	int projections = 10;
	int samples = 25;
	int horizon = 3;
	int depth = 5;
	real probability = (real) 4 / X.Rows();
	//real probability = 0.1;

	
	real bestIC50 = Y(0);
	for(int i=1;i<Y.Size();i++) {
		if(Y(i) > bestIC50)
			bestIC50 = Y(i);
	}
	real bestVal = 0;
	GaussianProcessTSBandit gpTS = GaussianProcessTSBandit(gp_scale_length,
		gp_scale_length_dic,gp_noise_var, gp_sig_var, 
		gp_threshold, delta);

	Vector predictions = Vector();
	Vector regret = regret.Null(iterations);
	Vector sregret = sregret.Null(iterations);
	bestVal = 0;
	gpTS.addObservation(X.getRow(42),Y(42));
	for(int i=0;i<iterations;i++) {
		int bestBandit = gpTS.predict(X, predictions, samples, horizon, probability);
		std::cout<<bestBandit<<",";
		gpTS.addObservation(X.getRow(bestBandit),Y(bestBandit));
		regret(i)=bestIC50 - Y(bestBandit);
		if(Y(bestBandit) > bestVal)
			bestVal = Y(bestBandit);
		sregret(i)=bestIC50 - bestVal;
	}	
	std::cout<<"avgCumSum=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<regret(i)<<" ";
	}
	std::cout<<"]";

	std::cout<<"avgSimRegret=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<sregret(i)<<" ";
	}
	std::cout<<"]";
	std::cout<<"\n -----\n"<<"GPTS regret: "<<regret.Sum()/iterations;
	std::cout<<"\n"<<"GPTS sregret: "<<sregret(iterations-1)<<"\n -----";
	std::cout<<"\n";

	GaussianProcessBandit gpUcb = GaussianProcessBandit(gp_scale_length,
		gp_scale_length_dic,gp_noise_var, gp_sig_var, 
		gp_threshold, delta);
	predictions = Vector();
	regret = regret.Null(iterations);
	sregret = sregret.Null(iterations);
	bestVal = 0;
	gpUcb.addObservation(X.getRow(42),Y(42));
	for(int i=0;i<iterations;i++) {
		int bestBandit = gpUcb.predict(X, predictions);
		gpUcb.addObservation(X.getRow(bestBandit),Y(bestBandit));
		regret(i)=bestIC50 - Y(bestBandit);
		if(Y(bestBandit) > bestVal)
			bestVal = Y(bestBandit);
		sregret(i)=bestIC50 - bestVal;
	}	
	std::cout<<"avgCumSum=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<regret(i)<<" ";
	}
	std::cout<<"]";

	std::cout<<"avgSimRegret=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<sregret(i)<<" ";
	}
	std::cout<<"]";
	std::cout<<"\n -----\n"<<"GPB regret: "<<regret.Sum()/iterations;
	std::cout<<"\n"<<"GPB sregret: "<<sregret(iterations-1)<<"\n -----";
	std::cout<<"\n";

	GaussianProcessRPBandit gpRPUcb = GaussianProcessRPBandit(gp_scale_length,
		gp_scale_length_dic,gp_noise_var, gp_sig_var, 
		gp_threshold, delta);
	predictions = Vector();
	regret = regret.Null(iterations);
	sregret = sregret.Null(iterations);
	bestVal = 0;
	gpRPUcb.addObservation(X.getRow(42),Y(42));
	for(int i=0;i<iterations;i++) {
		int bestBandit = gpRPUcb.predict(X, predictions,projections,depth);
		gpRPUcb.addObservation(X.getRow(bestBandit),Y(bestBandit));
		regret(i)=bestIC50 - Y(bestBandit);
		if(Y(bestBandit) > bestVal)
			bestVal = Y(bestBandit);
		sregret(i)=bestIC50 - bestVal;
	}
	std::cout<<"avgCumSum=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<regret(i)<<" ";
	}
	std::cout<<"]";

	std::cout<<"avgSimRegret=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<sregret(i)<<" ";
	}
	std::cout<<"]";
	std::cout<<"\n -----\n"<<"GPRPB regret: "<<regret.Sum()/iterations;
	std::cout<<"\n"<<"GPRPB sregret: "<<sregret(iterations-1)<<"\n -----";
	std::cout<<"\n";

	real R = 100;
	real epsilon = 0.5;

	LinearBanditTS lbTS = LinearBanditTS(R, delta, epsilon);
	predictions = Vector();
	regret = regret.Null(iterations);
	sregret = sregret.Null(iterations);
	bestVal = 0;
	lbTS.addObservation(X.getRow(42),Y(42));
	for(int i=0;i<iterations;i++) {
		int bestBandit = lbTS.predict(X, predictions);
		lbTS.addObservation(X.getRow(bestBandit),Y(bestBandit));
		regret(i)=bestIC50 - Y(bestBandit);
		if(Y(bestBandit) > bestVal)
			bestVal = Y(bestBandit);
		sregret(i)=bestIC50 - bestVal;
	}
	std::cout<<"avgCumSum=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<regret(i)<<" ";
	}
	std::cout<<"]";

	std::cout<<"avgSimRegret=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<sregret(i)<<" ";
	}
	std::cout<<"]";
	std::cout<<"\n -----\n"<<"LBTS regret: "<<regret.Sum()/iterations;
	std::cout<<"\n"<<"LBTS sregret: "<<sregret(iterations-1)<<"\n -----";
	std::cout<<"\n";

	regret = regret.Null(iterations);
	sregret = sregret.Null(iterations);
	bestVal = 0;
	for(int i=0;i<iterations;i++) {
		int index = urandom(0,Y.Size());
		regret(i)=bestIC50 - Y(index);
		if(Y(index) > bestVal)
			bestVal = Y(index);
		sregret(i)=bestIC50 - bestVal;
	}
	std::cout<<"avgCumSum=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<regret(i)<<" ";
	}
	std::cout<<"]";

	std::cout<<"avgSimRegret=[";
	for(int i=0;i<iterations;i++) {
		std::cout<<sregret(i)<<" ";
	}
	std::cout<<"]";
	std::cout<<"\n -----\n"<<"Random regret: "<<regret.Sum()/iterations;
	std::cout<<"\n"<<"Random sregret: "<<sregret(iterations-1)<<"\n -----";
	std::cout<<"\n";
}
