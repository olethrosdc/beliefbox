// -*- Mode: C++; -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

//#include "toy_uct_stopping_problem.h"

#include "UCB.h"
#include "Random.h"
#include "RandomNumberFile.h"
#include "ActionValueEstimate.h"
#include "DiscreteBanditPolicy.h"

/// A toy UCT stopping problem


//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int main (int argc, char** argv)
{
    real alpha = 1.0;
    real beta = 1.0;
    //int n_states = 1;

    int method = 0;
    if (argc > 1) {
        method = atoi(argv[1]);
    } else {
        fprintf(stderr, "Usage: bandit_ucb method gamma actions verbose experiments horizon\n");
    }

    real gamma = 0.99;
    if (argc > 2) {
        gamma = atof(argv[2]);
    }

    int n_actions = 2;
    if (argc > 3) {
        n_actions = atoi(argv[3]);
        assert(n_actions > 1);
    }

    int verbose = 0;
    if (argc > 4) {
        verbose = atoi(argv[4]);
    }

    int n_experiments = 1000;
    if (argc > 5) {
        n_experiments = atoi(argv[5]);
    }

    int horizon = (int) ceil(2.0/(1.0 - gamma));
    if (argc > 6) {
        horizon = atoi(argv[6]);
    }


    std::vector<real> rewards(horizon);
    for (int t=0; t<horizon; t++) {
        rewards[t] = 0.0;
    }
    
    real average_oracle_return = 0.0;
    real average_oracle_reward = 0.0;
    srand48(1228517343);
    //RandomNumberFile rng("./dat/r1e7.bin");
    for (int experiment=0; experiment<n_experiments; experiment++) {
        std::vector<real> Er(n_actions);

        for (int i=0; i<n_actions; i++) {
            Er[i] = drand48(); //rng.uniform();
        }

        DiscreteBanditPolicy* policy = NULL;
        BernoulliEstimate bernoulli_estimator(n_actions, 1, 1, 1);
        BernoulliEstimate biased_bernoulli_estimator(n_actions, 1, 1, 0);
        PointEstimate point_estimator(n_actions, 10.0 / (real) horizon, 1.0, 0.0);
        switch(method) {
        case 0:
            policy = new UCB1Policy(n_actions); 
            break;
        case 1:
            policy = new UCBgPolicy(n_actions, gamma);
            break;
        case 2:
            policy = new EpsilonGreedyPolicy(n_actions, 0.0, &bernoulli_estimator);
            break;
        case 3:
            policy = new EpsilonGreedyPolicy(n_actions, 0.0, &biased_bernoulli_estimator);
            break;
        case 4:
            policy = new EpsilonGreedyPolicy(n_actions, 0.0, &point_estimator);
            break;
        case 5:
            policy = new EpsilonGreedyPolicy(n_actions, 0.01, &point_estimator);
            break;
        case 6:
            policy = new EpsilonGreedyPolicy(n_actions, 0.1, &point_estimator);
            break;

        default:
            fprintf(stderr, "Unknown method %d\n", method);
            exit(-1);
        }

        for (int t=0; t<horizon; t++) {
            int action = policy->SelectAction();
            real reward = 0.0;
            
            if (urandom() < Er[action]) {
                reward = 1.0;
            } else {
                reward = 0;
            }
            
            rewards[t] += reward;
            policy->Observe(action, reward);
        }
        
        real oracle_return = 0.0;
        real max_reward = Max(Er);
        //printf ("%f\n", max_reward);
        oracle_return = max_reward*(1.0 - pow(gamma, (real) (horizon+1)))/(1.0 - gamma);
        
        average_oracle_return += oracle_return;
        average_oracle_reward += max_reward * ((real) horizon);

        delete policy;
    }


    real inv_exp = 1.0 / (real) n_experiments;
    real total_reward = 0.0;
    real discounted_reward = 0.0;
    real discount = 1.0;
    for (int t=0; t<horizon; t++) {
        //rewards[t] *= inv_exp;
        discounted_reward += discount * rewards[t];
        total_reward += rewards[t];
        discount *= gamma;
    }
    average_oracle_return *= inv_exp;
    average_oracle_reward *= inv_exp;
    discounted_reward *= inv_exp;
    total_reward *= inv_exp;



    printf("%d %f %d %f %f %f %f\n",
           method,
           gamma,
           n_actions,
           total_reward,
           discounted_reward,
           average_oracle_reward,
           average_oracle_return);
    return 0;
}



#endif
