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

#include "RSAPI.h"
#include "RandomNumberGenerator.h"

RolloutState::RolloutState(Environment<Vector, int>* environment_,
                           AbstractPolicy<Vector, int>* policy_,
						   Vector& start_state_,
                           real gamma_)
	: environment(environment_),
      policy(policy_),
	  start_state(start_state_),
      gamma(gamma_)
{
	V_U = 0;
    V_L = 0;
	V = 0;
    U = INF;
}

RolloutState::~RolloutState()
{
	for (uint i=0; i<rollouts.size(); ++i) {
		delete rollouts[i];
	}
}

void RolloutState::NewRollout(AbstractPolicy<Vector, int>* policy, int action)
{
	Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout
		= new Rollout<Vector, int, AbstractPolicy<Vector, int> >(start_state, action, policy, environment, gamma); //valgrind
	rollouts.push_back(rollout);
	
}	

/// Extend all rollouts by T
void RolloutState::ExtendAllRollouts(const int T)
{
	for (uint i=0; i<rollouts.size(); ++i) {
		rollouts[i]->Sample(T);
	}

}

/// Get the best empirical action in isolation
///
/// The error bounds used here are totally silly, though!
/// See OOP (Bubeck et al) for better bounds.
int RolloutState::BestEmpiricalAction()
{
    Vector Q_U(environment->getNActions()); // Upper bound on Q
    Vector Q_L(environment->getNActions()); // Lower bound on Q
    Vector N(environment->getNActions()); // number of times action was played
    real log_gamma = log(gamma);
	for (uint i=0; i<rollouts.size(); ++i) {
        Rollout<Vector, int, AbstractPolicy<Vector, int> >* rollout = rollouts[i];
        real error_bound = 0;
        if (rollout->running) {
            error_bound = exp(rollout->T * log_gamma);
        }
        int a = rollout->start_action;
        N(a)++;
        Q_U(a) += rollout->discounted_reward + error_bound;
        Q_L(a) += rollout->discounted_reward - error_bound;
    }
    Q_U = Q_U / N; // normalise
    Q_L = Q_L / N; // normalise
    
    policy->setState(start_state);
    int normal_action = policy->SelectAction();
    int optimistic_action = ArgMax(Q_U);
    int pessimistic_action = ArgMax(Q_L);
    if (Q_L(pessimistic_action) > Q_U(normal_action)) {
        return pessimistic_action;
    }  else {
        return normal_action;
    }
}

RSAPI::RSAPI(Environment<Vector, int>* environment_,
             RandomNumberGenerator* rng_,
             real gamma_)
    : environment(environment_),
      policy(NULL),
      rng(rng_),
      gamma(gamma_)
{
    // nothing to do here.
}


RSAPI::~RSAPI()
{
    for (uint i=0; i<states.size(); ++i) {
        delete states[i];
    }
}

void RSAPI::AddState(Vector& state)
{
    states.push_back(new RolloutState(environment, policy, state, gamma));
}

void RSAPI::SampleRandomly(const int T)
{
    RolloutState* state = states[rng->discrete_uniform(states.size())];
    state->ExtendAllRollouts(T);
}

void RSAPI::NewRandomRollouts(const int K, const int T)
{
    for (int k=0; k<K; ++k) {
        RolloutState* state = states[rng->discrete_uniform(states.size())];
        state->NewRollout(policy,
                          rng->discrete_uniform(environment->getNActions()));
        state->ExtendAllRollouts(T);
    }
}


/// Sample uniformly from all states
///
/// Take K rollouts from all state-action pairs, of length T
void RSAPI::SampleUniformly(const int K, const int T)
{                                                                       
    for (uint i=0; i<states.size(); ++i) {
        RolloutState* state = states[i];
        for (uint a=0; a<environment->getNActions(); ++a) {
            for (int k=0; k<K; ++k) {
                state->NewRollout(policy, a); // valgrind
            }
        }
        state->ExtendAllRollouts(T);
    }

}

/// Simply train the classifier
void RSAPI::TrainClassifier(Classifier<Vector, int, Vector>* classifier)
{
    for (uint i=0; i<states.size(); ++i) {
        int best_action = states[i]->BestEmpiricalAction();
        if (best_action >= 0) {
            //printf("# Action: %d state: ", best_action);
            //states[i]->start_state.print(stdout);
            classifier->Observe(states[i]->start_state, best_action);
        }
    }        
}

