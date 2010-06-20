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

#include "HQLearning.h"

HQLearning::HQLearning(int n_agents_,
					   int n_states_,
					   int n_actions_,
					   real gamma_,
					   real lambda_,
					   real alpha_,
					   real epsilon_, 
					   real beta_,
					   real initial_value_,
					   real baseline_)
    : n_agents(n_agents_),
	  n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      lambda(lambda_),
      alpha(alpha_),
	  epsilon(epsilon_),
	  beta(beta_),
      initial_value(initial_value_),
      baseline(baseline_),
      Q(n_states_, n_actions_),
      el(n_states_, n_actions_)
{
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q(s, a) = 0.0;//initial_value;
        }
    }
	
	sub_agent.resize(n_agents);
	for (int i=0; i<n_agents; ++i) {
		sub_agent[i] = new SubAgent(n_states, n_actions, gamma, lambda, alpha, epsilon, beta);
	}
    Reset();
}

void HQLearning::Reset()
{
	time_running = 0;
    state = -1;
    action = -1;
	previous_agent = 1;
	current_agent = 1;
	previous_observation = -1;
	current_observation = -1;
	terminating_observation = 0;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            el(s,a) = 0.0;
        }
    }
}

real HQLearning::Observe (int action, int next_state, real reward)
{
	return 0;
}

real HQLearning::Observe (real reward, int next_state, int next_action)
{
	return 0;
}

int HQLearning::Act(real reward, int next_state)
{
	previous_observation = current_observation;
	current_observation = next_state;
	time_running++;
	total_reward += reward;
	SubAgent* old_agent = sub_agent[current_agent];
	SubAgent* new_agent = old_agent;

	if (current_observation == terminating_observation) {
		previous_agent = current_agent;
		current_agent = (current_agent + 1)%n_agents;
		new_agent = sub_agent[current_agent];
		old_agent->HQ(terminating_observation)
			= (1 - alpha) * old_agent->HQ(terminating_observation)
			+ alpha * (total_reward + pow(gamma, time_running) * Max(new_agent->HQ));
		// Set up next agent.
		terminating_observation = new_agent->SelectTerminatingCondition();		
		time_running = 0;
		total_reward = 0;
		new_agent->Reset();
	}
 
    int next_action = new_agent->Act(reward, next_state);
    return next_action;
}
