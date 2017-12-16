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

#ifndef TREE_BRL_H
#define TREE_BRL_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "MDPModel.h"
#include "MultiMDPValueIteration.h"
#include "ValueIteration.h"
#include "PolicyEvaluation.h"
#include <vector>
#include <memory>

/// \ingroup ReinforcementLearning
/// @{
    
/** Direct model-based reinforcement learning using trees.
  
 */
class TreeBRL : public OnlineAlgorithm<int, int>
{
protected:
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real epsilon; ///< randomness
    int current_state; ///< current state
    int current_action; ///< current action
    MDPModel* belief; ///< pointer to the base MDP model
    RandomNumberGenerator* rng; ///< random number generator to draw samples from
    int horizon; ///< maximum number of samples to take
    int T; ///< time passed
    int size; ///< size of tree
public:
    class BeliefState
    {
    protected:
        TreeBRL& tree;
        MDPModel* belief;
        int state;
        int prev_action;
        real prev_reward;
        real probability;
        std::vector<BeliefState> children; ///< previous time
        BeliefState* prev; ///< previous belief state
        int t; ///< time
    public:
        /// This is used for the first belief state
        BeliefState(TreeBRL& tree_,
                    MDPModel* belief_,
                    int state_) : tree(tree_), belief(belief_), state(state_), probability(1), t(0)
        {
            tree.size++;
        }

        /// Use this to construct a subsequent belief state
        BeliefState(TreeBRL& tree_,
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
        void SparseExpandAllActions(int n_samples)
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
        void ExpandAllActions()
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

            
        
        /// Return the values using the mean MDP heuristic
        real CalculateValues()
        {
            Vector Q(tree.n_actions);
            //Q.Clear();
            //N.Clear();
            real V = prev_reward;
            if (t < tree.horizon) {
                for (uint i=0; i<children.size(); ++i) {
                    int a = children[i].prev_action;
                    Q(a) += children[i].probability * children[i].CalculateValues();
                }
                V += tree.gamma * Max(Q);
            } else {

                const DiscreteMDP* model = belief->getMeanMDP();
                ValueIteration VI(model, tree.gamma);
                VI.ComputeStateValuesStandard(1e-3);
                V += tree.gamma * VI.getValue(state);
            }
            //N.print(stdout);
							
            //printf("t: %d, r: %f, v: %f #MV\n", t, prev_reward, V);
            return V;
        }

        /// Return the values using an upper bound
        real CalculateUpperBoundValues(int n_samples)
        {
            Vector Q(tree.n_actions);
            real V = prev_reward;
            if (t < tree.horizon) {
                for (uint i=0; i<children.size(); ++i) {
                    int a = children[i].prev_action;
                    Q(a) += children[i].probability * children[i].CalculateUpperBoundValues(n_samples);
                }
                V += tree.gamma * Max(Q);
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
                V += tree.gamma * V_next / (real) n_samples;
            }
            //printf("t: %d, r: %f, v: %f, s:%d #UV\n", t, prev_reward, V, n_samples);
            return V;
        }

        /// Return the values using an upper bound
        real CalculateLowerBoundValues(int n_samples)
        {
            Vector Q(tree.n_actions);
            real V = prev_reward;
            if (t < tree.horizon) {
                for (uint i=0; i<children.size(); ++i) {
                    int a = children[i].prev_action;
                    Q(a) += children[i].probability * children[i].CalculateLowerBoundValues(n_samples);
                }
                V += tree.gamma * Max(Q);
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
                V += tree.gamma * V_next / (real) n_samples;
            }
            
            //printf("t: %d, r: %f, v: %f, s:%d #UV\n", t, prev_reward, V, n_samples);
            return V;
        }
        

        
    };
    TreeBRL(int n_states_, ///< number of states
            int n_actions_, ///< number of actions
            real gamma_, ///< discount factor
            MDPModel* belief_, ///< belief about the MDP
            RandomNumberGenerator* rng_, ///< the RNG
            int horizon_ = 1);
    virtual ~TreeBRL();
    virtual void Reset();
    virtual void Reset(int state);
    /// Full observation
    virtual real Observe (int state, int action, real reward, int next_state, int next_action);
    /// Partial observation 
    virtual real Observe (real reward, int next_state, int next_action);
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, int next_state);
    /** Set the rewards to Singular distributions.

        Since this is a Bayesian approach, we can simply set the belief about the reward in each state to be a singular distribution.
    */
    virtual void setFixedRewards(const Matrix& rewards)
    {
        belief->setFixedRewards(rewards);
#if 0
        logmsg("Setting reward matrix\n");
        rewards.print(stdout);
        belief->ShowModel();
#endif
    }

    virtual real getValue(int state, int action)
    {
        return 0;
    }

    /// Calculate a sparse belief tree where we take n_samples state
    /// samples and use n_TS MDP samples for the upper and lower bounds at
    /// the leaf nodes
    void CalculateSparseBeliefTree(int n_samples, int n_TS)
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
    void CalculateBeliefTree(int n_TS)
    {
        // Initialise the root belief state
        BeliefState belief_state(*this, belief, current_state);
        belief_state.ExpandAllActions();
        printf("%d %d %f %f %f\n", horizon, n_TS,
               belief_state.CalculateValues(),
               belief_state.CalculateLowerBoundValues(n_TS),
               belief_state.CalculateUpperBoundValues(n_TS));
    }
    
};


/// @}
#endif

