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
#include "Classifier.h"
#include "KDTree.h"
#include <vector>

class RandomNumberGenerator;

class RolloutState
{
public:
	Environment<Vector, int>* environment; ///< the environment
	AbstractPolicy<Vector, int>* policy; ///< the base policy to compare against
	Vector start_state; ///< The starting state
	real gamma;  ///< The value of gamma used
	std::vector<Rollout<Vector, int, AbstractPolicy<Vector, int> >* > rollouts; ///< The set of rollouts
	real V_U; ///< Upper bound on the value
	real V_L; ///< Lower bound on the value
	real V; ///< Best guess estimate
    real U; ///< The utility of sampling from that state
    RolloutState(Environment<Vector, int>* environment_,
                 AbstractPolicy<Vector, int>* policy_,
                 Vector& start_state_,
                 real gamma_);
	~RolloutState();
    Vector SampleFromPolicy();
	void NewRollout(AbstractPolicy<Vector, int>* policy, int action);
	void ExtendAllRollouts(const int T);
    Vector getRandomTerminalState();
    int BestEmpiricalAction();
    std::pair<Vector, bool> BestGroupAction();
    void Bootstrap(KDTree<RolloutState>& tree,
                   real L);
};


class RSAPI
{
protected:
    Environment<Vector, int>* environment; ///< The environment
    AbstractPolicy<Vector, int >* policy; ///< The current policy used
    RandomNumberGenerator* rng; ///< The random number generator
    real gamma; ///< Discount factor
public:
	std::vector<RolloutState*> states; ///< The set of states to sample from
    RSAPI(Environment<Vector, int>* environment_, RandomNumberGenerator* rng_, real gamma_);
    ~RSAPI();
    void AddState(Vector& state);
    Vector SampleStateFromPolicy() const;
    void setPolicy(AbstractPolicy<Vector, int>* policy_) 
    {
        assert(policy_);
        policy = policy_;
    }
    void SampleRandomly(const int T);
    void NewRandomRollouts(const int K, const int T);
    void SampleUniformly(const int K, const int T);
    int TrainClassifier(Classifier<Vector, int, Vector>* classifier);
    int GroupTrainClassifier(Classifier<Vector, int, Vector>* classifier);
    real LipschitzBound();
    void Bootstrap();
};


#endif
	
