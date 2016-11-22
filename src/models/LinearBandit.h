#ifndef LINEAR_BANDIT_H
#define LINEAR_BANDIT_H

#include <vector>
#include <random>
#include "real.h"
#include "Matrix.h"
#include "Vector.h"

class LinearBandit
{
private:
	Vector rewardData;
	Matrix examples;
	Vector rewardHistory;
	Vector predictionVector;
	Vector meanVector;
	std::vector<int> untested;
	std::vector<int> tested;

	Vector testedValues;
	Matrix testedData;

	real delta;
	int R;

	void play(int index);
	real calculateV(int T);

public:
	LinearBandit(Vector rewardData_, Matrix examples_, int R_, real delta_);
	void run(int times);
	~LinearBandit();
	void reset();
	Vector getPredictionVector();
	Vector getRewardHistory();
	Vector getMeanVector();
	std::vector<int> getTested();
};

#endif