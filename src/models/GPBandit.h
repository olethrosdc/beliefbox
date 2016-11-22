#ifndef GP_BANDIT_H
#define GP_BANDIT_H

#include <vector>
#include <random>
#include "GaussianProcess.h"
#include "real.h"
#include "Matrix.h"
#include "Vector.h"

class GPBandit
{
private:
	Vector rewardData;
	Matrix examples;
	Vector rewardHistory;
	GaussianProcess gp;
	Vector scale_length;
	Vector predictionVector;
	Vector meanVector;
	std::vector<int> untested;
	std::vector<int> tested;
	real delta;

	Vector testedValues;
	Matrix testedData;

	real betaFunction(real delta);
	void play(int index);

public:
	GPBandit(Vector rewardData_, Matrix examples_, Vector scale_length_, real noise_variance, real hyper_param, real delta_);
	void run(int times);
	~GPBandit();
	void reset();
	Vector getPredictionVector();
	Vector getRewardHistory();
	Vector getMeanVector();
	std::vector<int> getTested();
};

#endif