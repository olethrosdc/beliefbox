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

#ifndef HQ_LEARNING_H
#define HQ_LEARNING_H

#include "QLearning.h"
#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "Matrix.h"
#include "Vector.h"
#include "real.h"
#include "ExplorationPolicy.h"
#include "OnlineAlgorithm.h"
#include "Random.h"
#include <vector>

/** Hierarchical Q-learning.
	
	An online version of HQ-Learning, based on:
	
	Marco Wiering, Juergen Schmidhuber, 1997, "HQ-Learning", Adaptive Behaviour 6:2.
 */
class HQLearning : public OnlineAlgorithm<int,int>
{
public:
	class SubAgent
	{
	public:
		int n_states;
		int n_actions;
		real gamma;
		real lambda;
		real alpha;
		Vector HQ;
		real epsilon;
		real beta;
		QLearning* q_learning;
		VFExplorationPolicy* exploration_policy;
		SubAgent(int n_states_,
				 int n_actions_,
				 real gamma_,
				 real lambda_,
				 real alpha_,
				 real epsilon_,
				 real beta_)
			:  n_states(n_states_),
			   n_actions(n_actions_),
			   gamma(gamma_),
			   lambda(lambda_),
			   alpha(alpha_),
			   HQ(n_states),
			   epsilon(epsilon_),
			   beta(beta_)
		{
			exploration_policy = new MaxSoftmaxPolicy(n_actions, beta, epsilon);
			q_learning = new QLearning(n_states, n_actions, gamma, lambda, alpha, exploration_policy);
			for (int i=0; i<n_states; ++i) {
				HQ(i) = 0;
			}
		}
		~SubAgent()
		{
			delete q_learning;
			delete exploration_policy;
		}
		void Reset()
		{
			q_learning->Reset();
		}
		int Act(real reward, int state)
		{
			return q_learning->Act(reward, state);
		}
		int SelectTerminatingCondition()
		{
			if (urandom() < epsilon) {
				return ((int) floor(urandom() * n_states))%n_states;
			}
			return ArgMax(HQ);
		}
	};
protected:
	int time_running;
	real total_reward;
	int n_agents;
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real lambda; ///< eligibility trace decay rate
    real alpha; ///< learning rate 
	real epsilon; ///< probability with which best action won't be selected
	real beta; ///< nearly optimal values
    VFExplorationPolicy* exploration_policy; ///< exploration policy
    real initial_value; ///< initial value for Q values
    real baseline; ///< baseline reward
	std::vector<SubAgent*> sub_agent; ///< the sub-agent that will perform Q learning on just the observations.
    Matrix Q; ///< the matrix of Q-values
    Matrix el;

    int state; ///< current state
    int action; ///< current action
	int previous_agent;
	int current_agent;
	int previous_observation;
	int current_observation;
	int current_reward;
	int terminating_observation;
public:
    HQLearning(int n_agents,
			   int n_states_,
			   int n_actions_,
			   real gamma_,
			   real lambda_,
			   real alpha_,
			   real epsilon_,
			   real beta_,
			   real initial_value_= 0.0,
			   real baseline_ = 0.0);
    virtual ~HQLearning()
    {
		for (uint i=0; i<sub_agent.size(); ++i) {
			delete sub_agent[i];
		}
    }
    virtual void Reset();
    virtual real Observe (int action, int next_state, real reward);
    virtual real Observe (real reward, int next_state, int next_action);
    virtual int Act(real reward, int next_state);
    virtual real getValue (int s, int a)
    {
        return Q(state, action);
    }

    
};

#endif

