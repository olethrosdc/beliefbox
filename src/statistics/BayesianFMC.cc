#include "BayesianFMC.h"
BayesianFMC::BayesianFMC (int n_obs_, int n_actions_, int n_models_, real prior_)
	: n_obs(n_obs_),
	  n_actions(n_actions_),
	  n_models(n_models_),
	  prior(prior_)
{
	real p = 1.0 / (real) n_obs;
	for (int x=0; x<n_obs; ++x) {
		Pr_next[x] = p;
		for (int i=0; i<n_models; ++i) {
			Pr_obs(i, x) = p;
		}
	}
}

BayesianFMC::~BayesianFMC()
{
}

real BayesianFMC::Observe(int observation)
{
	real predicted_value = Pr_next[observation];
	
	real sum = 0;
	for (int i=0; i<n_models; ++i) {
		Pr[i] *= Pr_obs(i, observation);
		sum += Pr[i];
	}
	Pr /= sum;
	
	return predicted_value;
}
real BayesianFMC::Observe(int action, int observation);
real BayesianFMC::ObservationProbability (int action, int observation);
void BayesianFMC::Reset();
int BayesianFMC::predict(int a);
    
