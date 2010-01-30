#include "BayesianFMC.h"

BayesianFMC::BayesianFMC (int n_obs_, int n_actions_, int n_models_, real prior_)
	: n_obs(n_obs_),
	  n_actions(n_actions_),
	  n_models(n_models_),
	  prior(prior_),
      mc(n_models),
      Pr(n_models),
      logPr(n_models),
      Pr_next(n_obs)
{
	real p = 1.0 / (real) n_obs;
	for (int x=0; x<n_obs; ++x) {
		Pr_next[x] = p;
	}
    
    real sum = 0;
    for (int i=0; i<n_models; ++i) {
        mc[i] = new FactoredMarkovChain(n_actions, n_obs, i);
        Pr[i] = pow(prior, (real) i);
        sum += Pr[i];
    }
    for (int i=0; i<n_models; ++i) {
        Pr[i] /= sum;
        logPr[i] = log(Pr[i]);
    }

}

BayesianFMC::~BayesianFMC()
{
    for (int i=0; i<n_models; ++i) {
        delete mc[i];
    }
}

real BayesianFMC::Observe(int observation)
{
	real predicted_value = Pr_next[observation];
	
	real log_sum = LOG_ZERO;
	for (int i=0; i<n_models; ++i) {
		logPr[i] += log(mc[i]->Observe(observation));
		log_sum = logAdd(log_sum, logPr[i]);
	}
	logPr -= log_sum;
    Pr = exp(logPr);

	return predicted_value;
}

real BayesianFMC::Observe(int action, int observation)
{
    real predicted_value = 0;
	real log_sum = LOG_ZERO;
	for (int i=0; i<n_models; ++i) {
        real p = mc[i]->Observe(action, observation);
        real log_p = log(p);
        predicted_value += p * Pr[i];
		logPr[i] += log_p;
		log_sum = logAdd(log_sum, logPr[i]);
	}
	logPr -= log_sum;
    Pr = exp(logPr);

	return predicted_value;
}

real BayesianFMC::ObservationProbability (int action, int observation)
{
    real predicted_value = 0;
	for (int i=0; i<n_models; ++i) {
        real p = mc[i]->ObservationProbability(action, observation);
        predicted_value += p * Pr[i];
	}
    return predicted_value;
}

void BayesianFMC::Reset()
{
	for (int i=0; i<n_models; ++i) {
        mc[i]->Reset();
    }
}

int BayesianFMC::predict(int a)
{
    Serror("Not implemented\n");              
    return -1;
}
    
