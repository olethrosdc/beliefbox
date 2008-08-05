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

/// A toy UCT stopping problem


//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int main (int argc, char** argv)
{
    real alpha = 1.0;
    real beta = 1.0;
    int n_states = 2;
    int n_actions = 2;

    int method = 0;
    if (argc > 1) {
        method = atoi(argv[1]);
    }

    real gamma = 0.99;
    if (argc > 2) {
        gamma = atof(argv[2]);
    }

    int n_iter = 2;
    if (argc > 3) {
        n_iter = atoi(argv[3]);
    }

    int verbose = 0;
    if (argc > 4) {
        verbose = atoi(argv[4]);
    }

    int n_experiments = 1000;
    if (argc > 5) {
        n_experiments = atoi(argv[5]);
    }

    int horizon = 2.0/(1.0 - gamma);
    if (argc > 6) {
        horizon = atoi(argv[6]);
    }

    int value_iterations = n_iter + 1;
    if (argc > 7) {
        value_iterations = atoi(argv[7]);
    }

    

    std::vector<real> rewards(horizon);
    for (int t=0; t<horizon; t++) {
        rewards[t] = 0.0;
    }
    
    real average_oracle_return = 0.0;

    RandomNumberFile rng("./dat/r1e7.bin");
    for (int experiment=0; experiment<n_experiments; experiment++) {
        real actual_probability = rng.uniform(); //true_random(false); //urandom(0, 1);
        int state = 0;
        UCBPolicy* policy = NULL;
        if (method==0) {
            policy = new UCB1Policy(2); 
        } else {
            policy = new UCBgPolicy(2, gamma);
        }
        policy->SetTimes(1, INT_MAX);
        policy->SetReward(1, 0.5);
        for (int t=0; t<horizon; t++) {
            int action = policy->SelectAction();
            real reward = 0.0;
            
            if (state == 0 && action == 0) {
                if (urandom() < actual_probability) {
                    reward = 1.0;
                } else {
                    reward = -1.0;
                }
            }
            
            rewards[t] += reward;
            int next_state = 0;
            if (action == 1) {
                next_state = 1;
            }
            
            policy->Observe(action, (reward + 1.0) / 2.0);
            state = next_state;

            if (state==1) {
                break;
            }
        }
        
        real oracle_return = 0.0;
        if (actual_probability > 0.5) {
            oracle_return = (2.0 * actual_probability - 1.0)/(1.0 - gamma);
        }
        average_oracle_return += oracle_return;

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
    discounted_reward *= inv_exp;
    total_reward *= inv_exp;



    printf("%d %f %d %f %f %f\n",
           method,
           gamma,
           n_iter,
           total_reward,
           discounted_reward,
           average_oracle_return);
    return 0;
}



#endif
