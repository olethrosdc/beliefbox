// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SampleBasedRL.h"

SampleBasedRL::SampleBasedRL(int n_states_,
                             int n_actions_,
                             real gamma_,
                             real epsilon_,
                             MDPModel* model_,
                             RandomNumberGenerator* rng_,
                             int max_samples_,
                             bool use_upper_bound_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      epsilon(epsilon_),
      model(model_),
      tmpQ(n_actions),
      VU(n_states),
      VL(n_states),
      QU(n_states, n_actions),
      QL(n_states, n_actions),
      max_samples(max_samples_),
      rng(rng_),
      T(0),
      update_interval(1),
      next_update(0),
      use_upper_bound(use_upper_bound_),
      use_sampling_threshold(false),
      sampling_threshold(0.1),
      weights(max_samples)
{
    printf("# Starting Sample-Based-RL with %d samples, update interval %d\n",
           max_samples, update_interval);

    current_state = -1;

    real w_i = 1.0 / (real) max_samples;
    mdp_list.resize(max_samples);
    value_iteration.resize(max_samples);
    printf("# Generating mean MDP\n");
    //mdp_list[0] = model->getMeanMDP();
    for (int i=0; i<max_samples; ++i) {
        printf("# Generating sampled MDP\n");
        mdp_list[i] = model->generate();
        weights[i] = w_i;
        value_iteration[i] = new ValueIteration(mdp_list[i], gamma);
    }

    printf ("# Setting up MultiMPDValueIteration\n");
    multi_value_iteration = new MultiMDPValueIteration(weights, mdp_list, gamma);
    printf ("# Testing MultiMPDValueIteration\n");
    multi_value_iteration->ComputeStateActionValues(0,1);
}
SampleBasedRL::~SampleBasedRL()
{
    CalculateUpperBound(1e-6, 1e4);
    for (int s=0; s<n_states; ++s) {
        printf ("%f ", UpperBound(s));
    }
    printf(" # SampleBasedRL upper bound\n");
    CalculateLowerBound(1e-6,1e4);
    for (int s=0; s<n_states; ++s) {
        printf ("%f ", LowerBound(s));
    }
    printf(" # SampleBasedRL lower bound\n");
    for (int i=0; i<max_samples; ++i) {
        delete mdp_list[i];
        delete value_iteration[i];
    }
    delete multi_value_iteration;
}
void SampleBasedRL::Reset()
{
    current_state = -1;
    next_update = T;
    Resample();
    //model->Reset();
}
/// Full observation
real SampleBasedRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        model->AddTransition(state, action, reward, next_state);
    }
    return 0.0;
}
/// Partial observation 
real SampleBasedRL::Observe (real reward, int next_state, int next_action)
{
    if (current_state >= 0) {
        model->AddTransition(current_state, current_action, reward, next_state);
    }
    current_state = next_state;
    current_action = next_action;
    return 0.0;
}


void SampleBasedRL::Resample()
{
    for (int i=0; i<max_samples; ++i) {
        delete mdp_list[i];
        mdp_list[i] = model->generate();
#if 1
        logmsg("Generating MDP model %d\n", i);
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                printf ("%f ", mdp_list[i]->getExpectedReward(s, a));
            }
            printf("# state %d\n", s);
        }
#endif
    }
}


void SampleBasedRL::CalculateUpperBound(real accuracy, int iterations)
{
    for (int j=0; j<max_samples; ++j) {
        value_iteration[j]->setMDP(mdp_list[j]);
        value_iteration[j]->ComputeStateValuesStandard(accuracy, iterations);
    }
    
    real Z = 1.0 / (real) max_samples;
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            QU(s,a) = 0;
            for (int i=0; i<max_samples; i++) {
                QU(s, a) += value_iteration[i]->getValue(s, a);
            }
        }
    }
    QU *= Z;
    for (int s=0; s<n_states; ++s) {
#if 1
        VU(s) = QU(s, 0);
        for (int a=1; a<n_actions; ++a) {
            VU(s) = std::max<real>(VU(s), QU(s, a));
        }
#endif
    }
}

void SampleBasedRL::CalculateLowerBound(real accuracy, int iterations)
{
    multi_value_iteration->setMDPList(mdp_list);
    multi_value_iteration->ComputeStateValues(accuracy, iterations);
    
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            QL(s,a) = multi_value_iteration->getValue(s, a);
        }
    }
    for (int s=0; s<n_states; ++s) {
        VL(s) = QL(s, 0);
        for (int a=1; a<n_actions; ++a) {
            VL(s) = std::max<real>(VL(s), QL(s, a));
        }
    }
}


/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int SampleBasedRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);
    T++;

    // update the model
    if (current_state >= 0) {
        model->AddTransition(current_state, current_action, reward, next_state);
    }
    current_state = next_state;

    // Update MDPs
    //mdp_list[0] = model->getMeanMDP();
    // Do note waste much time generating MDPs
    
    bool do_update = false;
    if (use_sampling_threshold) {
        for (int i=0; i<max_samples; ++i) {
            real p = mdp_list[i]->getTransitionProbability(current_state, current_action, next_state);
            weights[i] *= p;
        }
        weights /= weights.Sum();
        if (Max(weights) > 1.0 - sampling_threshold) {
            //printf("Minimum: %f\n", Min(weights));
            //weights.print(stdout);
            do_update = true;
            for (int i=0; i<max_samples; ++i) {
                weights[i] = 1.0 / (real) max_samples;
            }
        }
    } else {
        if (T >= next_update) {    
            do_update = true;
        }
    }
    if (do_update) {    
        //printf("# update: %d\n", T);
        //model->ShowModel();
        update_interval += 1;//(int) (ceil)(1.01*(double) T);
        next_update = T + update_interval;
        Resample();
        if (use_upper_bound) {
            CalculateUpperBound(1e-6, 1e3);
        } else {
            CalculateLowerBound(1e-6, 1e3);
        }
    }

    // update values    
    if (use_upper_bound) {
        //CalculateUpperBound(0, 1);
        for (int i=0; i<n_actions; i++) {
            tmpQ[i] = UpperBound(next_state, i);
        }
    } else {
        //CalculateLowerBound(0, 1);
        for (int i=0; i<n_actions; i++) {
            tmpQ[i] = LowerBound(next_state, i);
        }

    }

    int next_action;
    real epsilon_t = epsilon / (1.0 + sqrt((real) T));

    // choose action
    if (urandom()<epsilon_t) {
        next_action = rng->discrete_uniform(n_actions);
    } else {
        next_action = ArgMax(tmpQ);
    }
    current_action = next_action;
    //printf("%f %d #epsilon\n", epsilon_t, action);
    
    return current_action;
}


