/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ROLLOUT_H
#define ROLLOUT_H

#include "Environment.h"
#include "AbstractPolicy.h"
#include "Random.h"

/** A rollout.
 */
template<typename S, typename A, typename P>
class Rollout
{
public:
	struct Observations{
		S state;
		A action;
		real reward;
		S next_state;
		bool endsim;
        Observations(const S& s, const A& a, real r, const S& s2, bool end)
            : state(s), action(a), reward(r), next_state(s2), endsim(end)
        {
        }
	};
	std::vector<Observations> esamples;
	std::vector<std::vector<Observations> > Samples;
	std::vector<real> total_rewards;
	std::vector<real> discounted_rewards;
//public:
	S start_state; ///< initial state
	A start_action; ///< initial action
	S end_state; ///< current final state
	P* policy; ///< policy to be used for continuation
	Environment<S, A>* environment; ///< environment model
	real gamma; ///< discount factor
	int T; ///< total steps used
	real total_reward; ///< total reward received
	real discounted_reward; ///< discounted reward received
	bool running; ///< whether the rollout needs to be sampled again
	bool sampling; ///< whether the rollout needs to samples the observations
	Rollout(const S& start_state_,
			const A& start_action_, 
			P* policy_, 
			Environment<S, A>* environment_,
			real gamma_) :
		start_state(start_state_),
		start_action(start_action_),
		policy(policy_),
		environment(environment_),
		gamma(gamma_)
	{
		end_state = start_state;
        T = 0;
		total_reward = 0;
		discounted_reward = 0;
		running = false;
		sampling = false;
	}
	Rollout(const S& start_state_,
			P* policy_,
			Environment<S, A>* environment_,
			real gamma_,
			bool sampling_) : 
		start_state(start_state_),
		policy(policy_),
		environment(environment_),
		gamma(gamma_),
		sampling(sampling_)
	{
		end_state = start_state;
		policy->setState(start_state);
		start_action = policy->SelectAction();
		T = 0;
		total_reward = 0;
		discounted_reward = 0;
		running = false;
	}
	/// Destructor
	~Rollout()
    {
    }
	
	void Act(const A& a, real discount_factor)
	{
		running = environment->Act(a);
		real reward = environment->getReward();
		S previous_state = end_state;
		//printf("A: %d, r: %f\n", a, reward);
		end_state = environment->getState();
		bool endsim = environment->getEndsim();

		if(sampling)
		{
			esamples.push_back(Observations(previous_state, a, reward, end_state, endsim));		
		}

		T++;
		total_reward += reward;
		discounted_reward += reward * discount_factor;
	}
	
	void Sample(const int period, bool use_start_state = true)
    {
		total_reward = 0;
		discounted_reward = 0;
		real discount_factor = 1;
        environment->Reset();
        if (use_start_state) {
            environment->setState(start_state);
        }
		running = true;
		int t = 0;
		//logmsg("Period, %d, start state: ", period); start_state.print(stdout);
		while (t < period || period < 0) {
            policy->setState(environment->getState());
			if (!running) {
//				printf("Period = %d , T = %d\n",period,t);
//				logmsg("Episode ended\n");
				break;
			}
			if (t==0) {
				if(sampling){
					esamples.clear();
					start_action = policy->SelectAction();
				}
				Act(start_action, discount_factor);
			} else {
				Act(policy->SelectAction(), discount_factor);
			}

			if (period >= 0) {
				discount_factor*= gamma;
			}
			++t;
			if (period < 0 && urandom() < 1.0 - gamma) {
				running = false;
            }
		}
        //logmsg("Stopping after %d steps\n", t);
		if(sampling){
			Samples.push_back(esamples);
			total_rewards.push_back(total_reward);
			discounted_rewards.push_back(discounted_reward);
		}
		//logmsg("Length: %d, R: %f, U: %f, State: ", T, total_reward, discounted_reward); end_state.print(stdout);
	}
    /// Always sample from the same starting state
	void Sampling(const int episodes, const int period)
	{
		for(int i = 0; i<episodes; ++i)
		{
			Sample(period, true);
		}
	}

    /// Sample uniformly within the environment bounds
	void UniformSampling(const int episodes, const int period)
	{
		for(int i = 0; i<episodes; ++i)
		{
			SetState(urandom(environment->StateUpperBound(), environment->StateLowerBound()));
			Sample(period, true);
		}
	}

    /// The starting state equal to the starting distribution
    void StartingDistributionSampling(const int episodes, const int period)
	{
		for(int i = 0; i<episodes; ++i)
		{
			SetState(urandom(environment->StateUpperBound(), environment->StateLowerBound()));
			Sample(period, true);
		}
	}


	void SetState(const S& start_state_)
	{
		start_state = start_state_;
	}
	S getState(const int& r, const int& s)
	{
		return Samples[r][s].state;
	}
	A getAction(const int& r, const int& s)
	{
		return Samples[r][s].action;
	}
	real getReward(const int& r, const int& s)
	{
		return Samples[r][s].reward;
	}
	S getNextState(const int& r, const int& s)
	{
		return Samples[r][s].next_state;
	}
	bool getEndsim(const int& r, const int& s)
	{
		return Samples[r][s].endsim;
	}
	int getNRollouts()
	{
		return Samples.size();
	}
	int getNSamples()
	{
		int n_samples = 0;
		for(int i = 0; i<getNRollouts(); ++i)
		{
			n_samples += Samples[i].size();
		}
		return n_samples;
	}
	int getNSamples(const int r)
	{
		return Samples[r].size();
	}
};


#endif
