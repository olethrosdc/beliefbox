#include "LinearBandit.h"
#include "MultivariateNormal.h"
#include "NormalDistribution.h"

LinearBandit::LinearBandit(Vector rewardData_, Matrix examples_, int R_, real delta_) 
{

	rewardData = rewardData_;
	examples = examples_;//.L2Norm();
	R = R_;
	delta = delta_;
	reset();

}
LinearBandit::~LinearBandit() {
}
void LinearBandit::run(int times) {
	assert(times<rewardData.Size());
	const int banditsToTry = 2;
	tested = std::vector<int>(0);
	std::cout << "Initializing LinearBandit and running " << times << " times.\n";
	std::cout << rewardData.Size() << " examples, " << examples.Columns() << " dimensions.\n";

	int d = examples.Columns();

	Vector meanHat = meanHat.Null(d);
	Vector f = f.Null(d);
	Matrix B = B.Null(d,d);
	for(int i=0;i<d;i++) {
		B(i,i) = 1;
	}

	for(int its=0;its<times;its++) {
		MultivariateNormal mvn = MultivariateNormal(meanHat,pow(calculateV(times),2)*B.Inverse());
		Vector meanSample = mvn.generate();

		real bestVal = Product(examples.getRow(untested[0]),meanSample);
		int bestIndex = 0;
		for(int i=1;i<untested.size();i++) {
			real tmp = Product(examples.getRow(untested[i]),meanSample);
			if(tmp>bestVal) {
				bestVal = tmp;
				bestIndex = i;
			}
		}
		play(bestIndex);
		real reward = rewardHistory(rewardHistory.Size()-1);
		Matrix tmp = tmp.Null(examples.Columns(),examples.Columns());
		MatrixProduct(examples.getRow(untested[bestIndex]),examples.getRow(untested[bestIndex]),tmp);
		B += tmp;
		f += examples.getRow(untested[bestIndex]) * reward;
		meanHat = B.Inverse() * f;
	}
}

real LinearBandit::calculateV(int T) {
	return R * sqrt(24 / (1/log((double)T)) * examples.Columns() * log(1/delta));
}
void LinearBandit::reset() {
	testedValues = Vector();
	rewardHistory = Vector();
	meanVector = Vector();
	testedData.Resize(1,examples.Columns());
	untested = std::vector<int>(rewardData.Size());
	for(int i=0;i<untested.size();i++) {
		untested[i] = i;
	}
	tested = std::vector<int>(0);
}
std::vector<int> LinearBandit::getTested() {
	return tested;
}

void LinearBandit::play(int index) {
	assert(index<untested.size());
	int realIndex = untested[index];
	tested.push_back(realIndex);
	// observe reward..
	Vector v = examples.getRow(realIndex);
	NormalDistribution n = NormalDistribution(0,1);
	real noise = n.generate();
	testedValues.AddElement(rewardData(realIndex)+noise);
	if(tested.size()==1)
	{
		testedData.setRow(0, v);
	}
	else {
		testedData = testedData.AddRow(v);
	}
	rewardHistory.AddElement(rewardData(realIndex)+noise);
}
Vector LinearBandit::getPredictionVector() {
	return predictionVector;
}
Vector LinearBandit::getRewardHistory() {
	return rewardHistory;
}
Vector LinearBandit::getMeanVector() {
	return meanVector;
}
