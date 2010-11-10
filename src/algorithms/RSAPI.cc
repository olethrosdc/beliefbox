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
#include "GeometricDistribution.h"

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
    Vector confidence_interval = pow(N, (real) -0.5);
    Q_U += confidence_interval;
    Q_L -= confidence_interval;
    V_U = Max(Q_U);
    V_L = Min(Q_L);
    V = 0.5 * (V_U + V_L);

    policy->setState(start_state);
    int normal_action = policy->SelectAction();
    int optimistic_action = ArgMax(Q_U);
    int pessimistic_action = ArgMax(Q_L);
    real gap = (Q_L(pessimistic_action) - Q_U(normal_action));
    if (gap > 0) {
        //printf ("# gap: %f, a: %d->%d, s: ", gap, normal_action, pessimistic_action);
        //start_state.print(stdout);

        return pessimistic_action;
    }  else {
        //printf ("# gap: %f, CI: ", gap); 
        //confidence_interval.print(stdout); 
        //printf("Q_U: "); Q_U.print(stdout); 
        //printf("Q_L: "); Q_L.print(stdout);
        //return normal_action;
        return normal_action;
    }
}


/** Sample from the current policy.

    This uses a gamma-discount sample from the current policy.
    The idea is the following.
    
    An infinite horizon problem with discount factor \f$\gamma\f$ is
    equivalent to a finite horizon problem with a random horizon,
    geometrically distributed, such that at any time \f$t > 0\f$, the
    probability that the horizon \f$T = t\f$ is equal to \f$1 - \gamma\f$.
    
    We can thus sample from that distribution easily.
 */
Vector RolloutState::SampleFromPolicy()
{

    GeometricDistribution horizon_distribution(1 - gamma);
    int horizon_sample = horizon_distribution.generate();
    policy->setState(start_state);
    int normal_action = policy->SelectAction();
    Rollout<Vector, int, AbstractPolicy<Vector, int> > rollout(start_state, normal_action, policy, environment, gamma);
    rollout.Sample(horizon_sample);
    return rollout.end_state;

}

/// Group best and worst actions together.
///
/// Give 0 probability to any completely dominated action.
/// Other actions get a probability proportional to the number
/// of other actions they dominate.
std::pair<Vector, bool> RolloutState::BestGroupAction()
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
    Vector confidence_interval = pow(N, (real) -0.5);
    Q_U += confidence_interval;
    Q_L -= confidence_interval;
    V_U = Max(Q_U);
    V_L = Min(Q_L);
    V = 0.5 * (V_U + V_L);

    policy->setState(start_state);

    std::pair<Vector, bool> retval;
    Vector& Pr = retval.first;
    Pr.Resize(environment->getNActions()); // Action probabilitie
    bool& domination = retval.second;
    domination = false;
    for (uint i=0; i<environment->getNActions(); ++i) {
        Pr(i) = 0.5;
        for (uint j=0; j<environment->getNActions(); ++j) {
            if (i==j) {
                continue;
            }
            if (Q_U(i) <= Q_L(j)) {
                Pr(i) = 0;
                domination = true;
            }
            if (Q_L(i) >= Q_U(j)) {
                Pr(i) += 1;
                domination = true;
            }
        }    
    }
    Pr /= Pr.Sum();
#if 0
    if (domination) {
        printf("# S: "); start_state.print(stdout);
        printf("# P: "); Pr.print(stdout);
    }
#endif
    return retval;
}

/** Bootstrap from values of neighbouring states.
    
    TODO This does not work at the moment.
 */

void RolloutState::Bootstrap(KDTree<RolloutState>& tree,
                             real L)
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
            Vector s_T = rollout->end_state;
            OrderedFixedList<KDNode> knn_list = tree.FindKNearestNeighbours(s_T, 3);
            std::list<std::pair<real, KDNode*> >::iterator knn_iter;
            for (knn_iter = knn_list.S.begin();
                 knn_iter != knn_list.S.end();
                 ++knn_iter) {
                KDNode* node = knn_iter->second;
            }

        }
        int a = rollout->start_action;
        N(a)++;
        Q_U(a) += rollout->discounted_reward + error_bound;
        Q_L(a) += rollout->discounted_reward - error_bound;
    }
    Q_U = Q_U / N; // normalise
    Q_L = Q_L / N; // normalise
    Vector confidence_interval = pow(N, (real) -0.5);
    Q_U += confidence_interval;
    Q_L -= confidence_interval;
    V_U = Max(Q_U);
    V_L = Min(Q_L);
    V = 0.5 * (V_U + V_L);


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

/// Sample a new state
Vector RSAPI::SampleStateFromPolicy() const
{
    int k = rng->discrete_uniform(states.size());
    return states[k]->SampleFromPolicy();
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
                state->NewRollout(policy, a); 
            }
        }
        state->ExtendAllRollouts(T);
    }

}

/// Simply train the classifier
int RSAPI::TrainClassifier(Classifier<Vector, int, Vector>* classifier)
{
    int n_improved_actions = 0;

    for (uint i=0; i<states.size(); ++i) {
        int best_action = states[i]->BestEmpiricalAction();
        if (best_action >= 0) {
            //printf("# Action: %d state: ", best_action);
            //states[i]->start_state.print(stdout)
            n_improved_actions++;
            classifier->Observe(states[i]->start_state, best_action);
        }
    }        
    return n_improved_actions;
}

real RSAPI::LipschitzBound()
{
    real L = 0;
    for (uint i=0; i<states.size(); ++i) {
        for (uint j=0; j<states.size(); ++j) {
            if (i==j) {
                continue;
            }
            real D_U = states[i]->V_U - states[j]->V_U;
            real D_L = states[i]->V_L - states[j]->V_L;
            real D = std::max(D_U, D_L);
            real s = EuclideanNorm(&states[i]->start_state,
                                   &states[j]->start_state);
            if (s > 0) {
                L = std::max(L, D / s);
            }
        }
    }
    return L;
}


void RSAPI::Bootstrap()
{
    /// Make a KNN tree.
    KDTree<RolloutState> tree(environment->getNStates());
    for (uint i=0; i<states.size(); ++i) {
        tree.AddVectorObject(states[i]->start_state, states[i]);
    }

    real L = LipschitzBound();
    for (uint i=0; i<states.size(); ++i) {
        states[i]->Bootstrap(tree, L);
    }
}


/// Simply train the classifier with action groups
int RSAPI::GroupTrainClassifier(Classifier<Vector, int, Vector>* classifier)
{
    int n_improved_actions = 0;
    for (uint i=0; i<states.size(); ++i) {
        std::pair<Vector, bool> retval = states[i]->BestGroupAction();
        if (retval.second) {
            n_improved_actions++;
            classifier->Observe(states[i]->start_state, retval.first);
        }                                
    }
    return n_improved_actions;
}

