#include "GPBandit.h"
#include "NormalDistribution.h"

GPBandit::GPBandit(Vector rewardData_, Matrix examples_, Vector scale_length_, real noise_variance, real hyper_param, real delta_) 
: gp(noise_variance, scale_length_, hyper_param)
{

	rewardData = rewardData_;
	examples = examples_;
	scale_length = scale_length_;
	delta = delta_;

	reset();

}
GPBandit::~GPBandit() {
}
void GPBandit::run(int times) {
	assert(times<rewardData.Size());
	const int banditsToTry = 1;
	tested = std::vector<int>(0);
	std::cout << "Initializing GPBandit and running " << times << " times.\n";
	std::cout << rewardData.Size() << " examples, " << examples.Columns() << " dimensions.\n";

	int n_untested = untested.size();
	for(int i=0;i<banditsToTry;i++) {
		int x = rand() % n_untested;
		play(x);
	}
	gp.Observe(testedData,testedValues);
	int its = 0;
	while(its+banditsToTry < times) {
		predictionVector = Vector(n_untested);
		meanVector = Vector(n_untested);
		real maxVal = 0;
		int bestIndex = -1;
		for(uint i=0;i<untested.size();i++) {
			real mean;
			real var;
			Vector z = examples.getRow(untested[i]);
			gp.Prediction(z,mean,var);
			real value = mean + sqrt(betaFunction(delta)) * var;
			if(value > maxVal) {
				maxVal = value;
				bestIndex = i;
			}
			predictionVector(i) = value;
			meanVector(i) = mean;
		}
		// play best bandit
		play(bestIndex);
		gp.Observe(testedData,testedValues);
		its++;
	}
	return;
}
void GPBandit::reset() {
	testedValues = Vector();
	rewardHistory = Vector();
	meanVector = Vector();
	testedData.Resize(1,examples.Columns());
	untested = std::vector<int>(rewardData.Size());
	for(uint i=0;i<untested.size();i++) {
		untested[i] = i;
	}
	tested = std::vector<int>(0);
}
std::vector<int> GPBandit::getTested() {
	return tested;
}

real GPBandit::betaFunction(real delta) {
	assert(delta>0 && delta<1);
	int dimensions = examples.Columns();
	int t = tested.size();
	return 2 * std::log(dimensions * (t * t) * (M_PI * M_PI) / (6 * delta));
}
void GPBandit::play(int index) {
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
Vector GPBandit::getPredictionVector() {
	return predictionVector;
}
Vector GPBandit::getRewardHistory() {
	return rewardHistory;
}
Vector GPBandit::getMeanVector() {
	return meanVector;
}
