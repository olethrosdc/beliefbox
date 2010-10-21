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

#ifndef RSAPI_H
#define RSAPI_H

#include "Environment.h"
#include "Rollout.h"
#include "Vector.h"
#include "AbstractPolicy.h"
#include <vector>

class RandomNumberGenerator;

class RolloutState
{
public:
	Environment<Vector, int>* environment;
	Vector start_state; ///< The starting state
	real gamma;  ///< The value of gamma used
	std::vector<Rollout<Vector, int, AbstractPolicy<Vector, int> >* > rollouts; ///< The set of rollouts
	real V_U; ///< Upper bound on the value
	real V_L; ///< Lower bound on the value
	real V; ///< Best guess estimate
	RolloutState(Environment<Vector, int>* environment_,
				 Vector& start_state_);
	~RolloutState();
	void NewRollout(AbstractPolicy<Vector, int>* policy, int action);
	void Sample(const int T);
    Vector getRandomTerminalState();
};


class RSAPI
{
protected:
    Environment<Vector, int>* environment;
    AbstractPolicy<Vector, int >* policy;
    RandomNumberGenerator* rng;
public:
	std::vector<RolloutState*> states;
    RSAPI(Environment<Vector, int>* environment_, RandomNumberGenerator* rng_);
    ~RSAPI();
    void AddState(Vector& state);
    void setPolicy(AbstractPolicy<Vector, int>* policy_) 
    {
        assert(policy_);
        policy = policy_;
    }
    void SampleRandomly(const int T);
    void NewRandomRollouts(const int K, const int T);

};


#endif
	
