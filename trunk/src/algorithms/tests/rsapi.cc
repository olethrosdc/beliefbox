/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "RSAPI.h"
#include "Rollout.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "RandomPolicy.h"
#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
#include "KNNClassifier.h"
#include "ClassifierPolicy.h"

struct PerformanceStatistics
{
    real total_reward;
    real discounted_reward;
    real run_time;
    PerformanceStatistics()
        : total_reward(0),
          discounted_reward(0),
          run_time(0)
    {}
};

struct AveragePerformanceStatistics : PerformanceStatistics
{
    int N;
    AveragePerformanceStatistics()
        : N(0)
    {
    }
    void Observe(PerformanceStatistics& x)
    {
        N++;
        real rN = (real) N;
        real irN = 1 / rN;
        real gamma = (rN - 1) * irN;

        total_reward = gamma * total_reward + (1 - gamma) * x.total_reward;
        discounted_reward = gamma * discounted_reward + (1 - gamma) * x.discounted_reward;
        run_time = gamma * run_time + (1 - gamma) * x.run_time;

    }
    void Show()
    {
        printf ("# UR: %f, DR: %f, T: %f\n", 
                total_reward,
                discounted_reward,
                run_time);
    }
};
    

PerformanceStatistics Evaluate(Environment<Vector, int>* environment,
                               AbstractPolicy<Vector, int>* policy,
                               real gamma, int T);

int main(void)
{
	MersenneTwisterRNG rng;

	// Create a new environment
	Environment<Vector, int>* environment;

	//environment = new MountainCar();
	environment = new Pendulum();
    
	// Place holder for the policy
	AbstractPolicy<Vector, int>* policy;
	
	// Start with a random policy!
	policy = new RandomPolicy(environment->getNActions(), &rng);
	


    int n_states = 100;
    
    int state_dimension = environment->getNStates();
    Vector S_L = environment->StateLowerBound();
    Vector S_U = environment->StateUpperBound();
    
    printf("# State dimension: %d\n", state_dimension);
    printf("# S_L: "); S_L.print(stdout);
    printf("# S_U: "); S_U.print(stdout);
    KNNClassifier* classifier = NULL;
    real gamma = 0.99;

    std::vector<Vector> state_vector(n_states);
    for (int k=0; k<n_states; ++k) {
        Vector& state = state_vector[k];
        state.Resize(S_L.Size());
        for (int i=0; i<S_L.Size(); ++i) {
            state(i) = rng.uniform(S_L(i), S_U(i));
        }
    }
    int n_iter=10;
    for (int iter=0; iter<n_iter; ++iter) {
        AveragePerformanceStatistics statistics;
        for (int i=0; i<100; ++i) {
            PerformanceStatistics run_statistics = Evaluate(environment,
                                                            policy,
                                                            gamma,
                                                            1000);
            statistics.Observe(run_statistics);
        }
        statistics.Show();
        RSAPI rsapi(environment, &rng, gamma);
        rsapi.setPolicy(policy);
        for (int k=0; k<n_states; ++k) {
            rsapi.AddState(state_vector[k]);
        }
        
        rsapi.SampleUniformly(100,1000);

        KNNClassifier* new_classifier = new KNNClassifier(state_dimension, environment->getNActions(), 1);
        rsapi.TrainClassifier(new_classifier);
        delete policy;
        policy = new ClassifierPolicy(new_classifier);
        if (classifier) {
            delete classifier;
        }        
        classifier = new_classifier;
    }
    delete classifier;
    delete policy;
	delete environment;
}


PerformanceStatistics Evaluate(Environment<Vector, int>* environment,
                               AbstractPolicy<Vector, int>* policy,
                               real gamma, int T)
{
    PerformanceStatistics statistics;
    statistics.total_reward = 0;
    statistics.discounted_reward = 0;
    statistics.run_time = 0;

    real discount = 1;
    environment->Reset();
    for (int t=0; t < T; ++t, ++statistics.run_time) {
        Vector state = environment->getState();
        real reward = environment->getReward();
        statistics.total_reward += reward;
        statistics.discounted_reward += reward * discount;
        discount *= gamma;
        policy->setState(state);
        int action = policy->SelectAction();
        bool action_ok = environment->Act(action);
        if (!action_ok) {
            break;
        }
    }
    return statistics;
}
#endif
