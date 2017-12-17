// -*- Mode: c++ -*-
// copyright (c) 2017 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "TreeBRL.h"

TreeBRL::TreeBRL(int n_states_,
                 int n_actions_,
                 real gamma_,
                 MDPModel* belief_,
                 RandomNumberGenerator* rng_,
                 int horizon_)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      belief(belief_),
      rng(rng_),
      horizon(horizon_),
      T(0),
      size(0),
      Qs(n_actions)
{
  
    //printf("# Starting Tree-Bayes-RL with %d states, %d actions, %d horizon\n", n_states, n_actions, horizon);

    current_state = -1;

}

TreeBRL::~TreeBRL()
{
    //printf(" # destroying tree of size %d\n", size);
}

void TreeBRL::Reset()
{
    current_state = -1;
}

void TreeBRL::Reset(int state)
{
    current_state = state;
}

/// Full observation
real TreeBRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        belief->AddTransition(state, action, reward, next_state);
    }
    current_state = next_state;
    current_action = next_action;
    return 0.0;
}
/// Partial observation 
real TreeBRL::Observe (real reward, int next_state, int next_action)
{
    if (current_state >= 0) {
        belief->AddTransition(current_state, current_action, reward, next_state);
    }
    current_state = next_state;
    current_action = next_action;
    return 0.0;
}


/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int TreeBRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);
    T++;
    int n_MDP_leaf_samples = 2;
    CalculateBeliefTree(n_MDP_leaf_samples);
    int next_action = ArgMax(Qs);
    Observe(reward, next_state, next_action);
    current_action = next_action;
    current_state = next_state;
    return current_action;
}


    /// Calculate a sparse belief tree where we take n_samples state
    /// samples and use n_TS MDP samples for the upper and lower bounds at
    /// the leaf nodes

void TreeBRL::CalculateSparseBeliefTree(int n_samples, int n_TS)
{
    // Initialise the root belief state
    BeliefState belief_state(*this, belief, current_state);
    belief_state.SparseExpandAllActions(n_samples);
    printf("%d %d %f %f %f\n", horizon, n_TS,
           belief_state.CalculateValues(),
           belief_state.CalculateLowerBoundValues(n_TS),
           belief_state.CalculateUpperBoundValues(n_TS));
}

/// Calculate a sparse belief tree where we take n_TS MDP samples
/// for the upper and lower bounds at the leaf nodes
void TreeBRL::CalculateBeliefTree(int n_TS)
{
    // Initialise the root belief state
    BeliefState belief_state(*this, belief, current_state);
    belief_state.ExpandAllActions();
    printf("%d %d %f %f %f\n", horizon, n_TS,
           belief_state.CalculateValues(),
           belief_state.CalculateLowerBoundValues(n_TS),
           belief_state.CalculateUpperBoundValues(n_TS));
}

//------------- Belief states ----------------//

TreeBRL::BeliefState::BeliefState(TreeBRL& tree_,
                                  MDPModel* belief_,
                                  int state_) : tree(tree_), belief(belief_), state(state_), probability(1), t(0)
{
    tree.size++;
}

/// Use this to construct a subsequent belief state
TreeBRL::BeliefState::BeliefState(TreeBRL& tree_,
                                  MDPModel* belief_,
                                  int prev_state_,
                                  int prev_action_,
                                  int state_,
                                  real r,
                                  real p,
                                  BeliefState* prev_) : tree(tree_), belief(belief_), state(state_), prev_action(prev_action_), prev_reward(r), probability(p), prev(prev_), t(prev_->t + 1)
{
    belief->AddTransition(prev_state_, prev_action, r, state);
    tree.size++;
}
/// Generate transitions from the current state for all
/// actions. Do this recursively, using the marginal
/// distribution, but using sparse sampling.
void TreeBRL::BeliefState::SparseExpandAllActions(int n_samples)
{
    if (t >= tree.horizon) {
        return;
    }
    real p = 1 / (real) n_samples;
    for (int k=0; k<n_samples; ++k) {
        for (int a=0; a<tree.n_actions; ++a) {
            int next_state = belief->GenerateTransition(state, a);
            real reward = belief->GenerateReward(state, a);
            // Generate the new belief state and put it in the tree
            children.push_back(BeliefState(tree, belief, state, a, next_state, reward, p, this));
        }
    }

    for (uint i=0; i<children.size(); ++i) {
        children[i].SparseExpandAllActions(n_samples);
    }
}
/// Generate transitions from the current state for all
/// actions. Do this recursively, using the marginal distribution. Only expand the children when we're under the horizon.
void TreeBRL::BeliefState::ExpandAllActions()
{
    if (t >= tree.horizon) {
        return;
    }
    for (int a=0; a<tree.n_actions; ++a) {
        for (int next_state=0;
             next_state<tree.n_states;
             ++next_state) {
            real p = belief->getTransitionProbability(state, a, next_state);
            real reward = belief->GenerateReward(state, a);
            children.push_back(BeliefState(tree, belief, state, a, next_state, reward, p, this));
        }
    }
            
    for (uint i=0; i<children.size(); ++i) {
        children[i].ExpandAllActions();
    }
}

            
        
/// Return the values using the mean MDP heuristic.
///
/// \f$V_t(s) = r(s'',a'',s) + \gamma \max_a Q_{t+1} (s,a)\f$
real TreeBRL::BeliefState::CalculateValues()
{
    Vector Q(tree.n_actions);
    //Q.Clear();
    //N.Clear();
    real V = prev_reward;
    real discount = tree.gamma;
    if (t < tree.horizon) {
        for (uint i=0; i<children.size(); ++i) {
            int a = children[i].prev_action;
            Q(a) += children[i].probability * children[i].CalculateValues();
        }
        V += discount * Max(Q);
    } else {
        const DiscreteMDP* model = belief->getMeanMDP();
        ValueIteration VI(model, tree.gamma);
        VI.ComputeStateValuesStandard(1e-3);
        V += VI.getValue(state);
    }
    //N.print(stdout);

    tree.Qs = Q;
    //printf("t: %d, r: %f, v: %f #MV\n", t, prev_reward, V);
    return V;
}

/// Return the values using an upper bound
real TreeBRL::BeliefState::CalculateUpperBoundValues(int n_samples)
{
    Vector Q(tree.n_actions);
    real V = prev_reward;
    real discount = tree.gamma;
    if (t < tree.horizon) {
        for (uint i=0; i<children.size(); ++i) {
            int a = children[i].prev_action;
            Q(a) += children[i].probability * children[i].CalculateUpperBoundValues(n_samples);
        }
        V += discount * Max(Q);
    } else {
        //printf("t >= T : %d\n", n_samples);
        real V_next = 0;
        for (int i=0; i<n_samples; ++i) {
            DiscreteMDP* model = belief->generate();
            ValueIteration VI(model, tree.gamma);
            VI.ComputeStateValuesStandard(1e-3);
            V_next += VI.getValue(state);
            delete model;
        }
        V +=  V_next / (real) n_samples;
    }
    //printf("t: %d, r: %f, v: %f, s:%d #UV\n", t, prev_reward, V, n_samples);
    tree.Qs = Q;
    return V;
}

/// Return the values using an upper bound
real TreeBRL::BeliefState::CalculateLowerBoundValues(int n_samples)
{
    Vector Q(tree.n_actions);
    real V = prev_reward;
    real discount = tree.gamma;
    if (t < tree.horizon) {
        for (uint i=0; i<children.size(); ++i) {
            int a = children[i].prev_action;
            Q(a) += children[i].probability * children[i].CalculateLowerBoundValues(n_samples);
        }
        V += discount * Max(Q);
    } else {
        //printf("t >= T : %d\n", n_samples);
        const DiscreteMDP* model = belief->getMeanMDP();
        ValueIteration VI(model, tree.gamma);
        VI.ComputeStateValuesStandard(1e-3);
        FixedDiscretePolicy* policy = VI.getPolicy();
        real V_next = 0;
        for (int i=0; i<n_samples; ++i) {
            DiscreteMDP* model = belief->generate();
            PolicyEvaluation PI(policy, model, tree.gamma);
            PI.ComputeStateValues(1e-3);
            V_next += PI.getValue(state);
            delete model;
        }
        V += V_next / (real) n_samples;
    }
            
    //printf("t: %d, r: %f, v: %f, s:%d #UV\n", t, prev_reward, V, n_samples);
    tree.Qs = Q;
    return V;
}
        

        
