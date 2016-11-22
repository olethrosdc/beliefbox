#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "ReadData.h"
#include "Util.h"
#include "LinearBandit.h"
#include "MultivariateNormal.h"
#include "Dirichlet.h"

int main(int argc, char* argv[]) {
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
	//
	/*
	// NO COMPRESS
	for(int i=0;i<compressed.size();i++) {
		compressed[i].resize(c.indexMap.size());
		for(int j=0;j<descriptorIndex[i].size();j++) {
			int realIndex = descriptorIndex[i][j] - 1;

			compressed[i][descriptorIndex[i][j]-1] = descriptorCount[i][j];
		}
	}*/
	//

	// recommended by Lars
	//real scale_length = 2 * std::pow(10.0,-4);
	real scale_length = 2 * std::pow(10.0,2);
	real noise_var = 1.0;
	real hyper_param = 1.0;

	int iterations = 5;
	int samples = 1000;

	Vector rewardData_ = Vector(rewardData);
	Matrix compressed_ = Matrix(compressed.size(),compressed[0].size());
	for(int i=0;i<compressed.size();i++) {
		compressed_.setRow(i,Vector(compressed[i]));
	}
	//int b25 = samples * 0.25; //use GP to generate data
	int b25 = rewardData_.Size() * 0.25;
	Vector avgCumSum = Vector(b25);
	Vector cumSum = Vector(b25);
	Vector simRegret = Vector(b25);
	Vector avgSimRegret = Vector(b25);

	for(int t = 0;t<iterations;t++) {
		Matrix newPoints_ = Matrix(samples,compressed[0].size());
		Vector newRewards_ = Vector(samples);
		//generateData(compressed_,rewardData_,newPoints_,newRewards_,samples); // use GP to generate data

		newPoints_ = compressed_;
		newRewards_ = rewardData_;

		Vector v = v.Null(newPoints_.Columns())+scale_length;
		GPBandit gp = GPBandit(newRewards_, newPoints_, v, noise_var, hyper_param, 0.5);
		//LinearBandit gp = LinearBandit(newRewards_, newPoints_, 1000, 0.01);

		std::vector<real> newRewards(newRewards_.Size());
		for(int i = 0;i<newRewards.size();i++) {
			newRewards[i] = newRewards_(i);
		}

		gp.run(b25);
		std::vector<int> tested = gp.getTested();

		//Vector optVector = Vector(findNMaximum(newRewards,b25));
		Vector rewardHistory = gp.getRewardHistory();
		Vector meanVector = meanVector.Null(newRewards_.Size());
		for(int j=0;j<rewardHistory.Size();j++) {
			meanVector(tested[j]) += rewardHistory(j);
		}
		for(int j=0;j<newRewards_.Size();j++) {
			if(find(tested.begin(),tested.end(),j)!=tested.end())
				meanVector(j) /= std::count(tested.begin(),tested.end(),j);
		}
		real maxValue = meanVector(0);
		for(int j=1;j<meanVector.Size();j++) {
			if(meanVector(j) > maxValue) 
				maxValue = meanVector(j);

		}
		//Vector optVector = Vector(findNMaximum(newRewards,b25));
		//Vector optVector = optVector.Null(b25) + maxValue;

		real opt = Vector(findNMaximum(newRewards,1))(0);
		Vector optVector = optVector.Null(b25) + opt;
		real last = 0;
		real bestSoFar =0;
		for(int j=0;j<b25;j++) {
			if(newRewards_(tested[j]) > bestSoFar) {
				bestSoFar = newRewards_(tested[j]);
			}
			cumSum(j) = last + (optVector(j) - newRewards_(tested[j]));
			last = cumSum(j);
			cumSum(j) /= j+1;
			avgSimRegret(j) = optVector(0) - bestSoFar;
		}
		avgCumSum += cumSum;
		simRegret += avgSimRegret;
		gp.reset();

	}
	std::cout<<"avgCumSum=[";
	avgCumSum /= iterations;
	for(int i=0;i<b25;i++) {
		std::cout<<avgCumSum(i)<<" ";
	}
	std::cout<<"]";

	std::cout<<"avgSimRegret=[";
	simRegret /= iterations;
	for(int i=0;i<b25;i++) {
		std::cout<<simRegret(i)<<" ";
	}
	std::cout<<"]";
	return 0;
}
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
void generateData(Matrix oldPoints, Vector oldRewards, Matrix &newPoints, Vector &newRewards, int numOfPoints) {
	real scale_length = 2 * std::pow(10.0,2);
	real noise_var = 1.0;
	real hyper_param = 1.0;
	Vector v = v.Null(oldPoints.Columns())+scale_length;

	GaussianProcess gp(noise_var,v,hyper_param);
	DirichletDistribution dir = DirichletDistribution(oldPoints.Rows(), 1.0);

	// Learn GP 
	gp.Observe(oldPoints, oldRewards);
	for(int i=0;i<numOfPoints;i++) {
		Vector sample = dir.generate();
		Vector point = point.Null(oldPoints.Columns());
		for(int j=0;j<oldPoints.Rows();j++) {
			real s = sample(j);
			point+= oldPoints.getRow(j) * s;
		}
		real mean;
		real var;
		gp.Prediction(point,mean,var);
		newPoints.setRow(i,point);
		newRewards(i) = mean;
	}
}
