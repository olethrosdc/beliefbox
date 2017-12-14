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
#include <vector>

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
public:
    class BeliefState
    {
    protected:
        TreeBRL& tree;
        MDPModel* belief;
        int state;
        int prev_action;
        int prev_reward;
        std::vector<BeliefState> children; ///< previous time
        BeliefState* prev; ///< previous belief state
        int t; ///< time
    public:
        /// This is used for the first belief state
        BeliefState(TreeBRL& tree_,
                    MDPModel* belief_,
                    int state_) : tree(tree_), belief(belief_), state(state_), t(0)
        {}

        /// Use this to construct a subsequent belief state
        BeliefState(TreeBRL& tree_,
                    MDPModel* belief_,
                    int prev_state_,
                    int prev_action_,
                    int state_,
                    real r,
                    BeliefState* prev_) : tree(tree_), belief(belief_), state(state_), prev_action(prev_action_), prev_reward(r), prev(prev_), t(prev_->t + 1)
        {
            belief->AddTransition(prev_state_, prev_action, r, state);
        }
        /// Generate transitions from the current state for all
        /// actions. Do this recursively, using the marginal distribution
        void ExpandAllActions()
        {
            for (int a=0; a<tree.n_actions; ++a) {
                int next_state = belief->GenerateTransition(state, a);
                real reward = belief->GenerateReward(state, a);
                children.push_back(BeliefState(tree, belief, state, a, next_state, reward, this));
            }
            if (t < tree.horizon) {
                for (uint i=0; i<children.size(); ++i) {
                    children[i].ExpandAllActions();
                }
            }
        }

        /// place holder
        real CalculateValues()
        {
            real y = 0;
            if (t < tree.horizon) {
                for (uint i=0; i<children.size(); ++i) {
                    y += children[i].CalculateValues();
                }
            }
            return y;
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

    void CalculateBeliefTree()
    {
        // Initialise the root belief state
        BeliefState belief_state(*this, belief, current_state);
        belief_state.ExpandAllActions();
        belief_state.CalculateValues();
    }
};


/// @}
#endif

