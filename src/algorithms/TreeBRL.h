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
public:
	enum LeafNodeValue {
		NONE = 0x0, V_MIN, V_MAX, V_MEAN, V_UTS, V_LTS
	};
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
	int state_samples; ///< number of states to sample when branching
	int policy_samples; ///< number of actions (or policies) to sample when branching
	int reward_samples; ///< number of reward values to sample when branching
    int T; ///< time passed
    int size; ///< size of tree
    Vector Qs; ///< caching the value of the actions for the current state
	LeafNodeValue leaf_node_expansion; ///< how to expand the leaf node
public:
    class BeliefState
    {
    protected:
        TreeBRL& tree; ///< link to the base tree
		MDPModel* belief; ///< current belief
        int state; ///< current state
        int prev_action; ///< action taken to arrive here
        real prev_reward; ///< reward received to arrive here
        real probability; ///< probability of arriving here given previous state and action
        std::vector<BeliefState*> children; ///< next belief states
        BeliefState* prev; ///< previous belief state
        int current_step; ///< time
    public:
        BeliefState(TreeBRL& tree_,
                    const MDPModel* belief_,
                    int state_);
        BeliefState(TreeBRL& tree_,
                    const MDPModel* belief_,
                    int prev_state_,
                    int prev_action_,
                    int state_,
                    real r,
                    real p,
                    BeliefState* prev_);
		~BeliefState();
        // methods for building the tree
        void ExpandAllActions();
        void SparseExpandAllActions(int n_samples);
        // methods for calculating action values in the tree
        real CalculateValues(LeafNodeValue leaf_node);
		real MeanMDPValue();
        real UTSValue();
        real LTSValue();
		void print() const;
        // methods for adaptively building the tree while calculating values (TODO)
        // real StochasticBranchAndBound(int n_samples);
    };
    TreeBRL(int n_states_, ///< number of states
            int n_actions_, ///< number of actions
            real gamma_, ///< discount factor
            MDPModel* belief_, ///< belief about the MDP
            RandomNumberGenerator* rng_, ///< the RNG
            int horizon_ = 1,
			int state_samples_ = 2,
			int policy_samples_ = 2,
			int reward_samples_ = 2,
			LeafNodeValue leaf_node_expansion = NONE);
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

        Since this is a Bayesian approach, we can simply set the belief about the reward in each state to be a singular distribution, if we want. This would correspond to us having a fixed, unshakeable belief about them.
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
		Serror("Not implemented in this context\n");
        return 0;
    }

    TreeBRL::BeliefState CalculateSparseBeliefTree(int n_samples, int n_TS);
	TreeBRL::BeliefState CalculateBeliefTree();

    
};


/// @}
#endif

