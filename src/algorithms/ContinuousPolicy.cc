// -*- Mode: c++ -*-
// (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
// (c) 2021 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "ContinuousPolicy.h"
#include "Random.h"
#include <cmath>

///Create an uniform policy with a given number of basis functions.
FixedContinuousPolicy::FixedContinuousPolicy(int n_dimension, int n_actions_, BasisSet<Vector, int>& bfs_)
	: ContinuousPolicy(n_dimension, n_actions_, bfs_ ),
      epsilon_greedy(false),
      epsilon(0.0)
{
	logmsg("Initialising epsilon-greedy continuous policy");
	Vector w(n_actions*(bfs_.size() + 1));
	weights = w;
	p.Resize(n_actions);
}

FixedContinuousPolicy::FixedContinuousPolicy(int n_dimension_, int n_actions_, BasisSet<Vector, int>& bfs_, const Vector& weights_)
	: ContinuousPolicy(n_dimension_, n_actions_, bfs_, weights_),
      epsilon_greedy(false),
      epsilon(0.0)
{
	logmsg("Initialising epsilon-greedy continuous policy");
	p.Resize(n_actions);
}

/// Destructor. Does nothing at the moment.
FixedContinuousPolicy::~FixedContinuousPolicy()
{
}

int FixedContinuousPolicy::SelectAction()
{
	StatePolicy();
    if (epsilon_greedy && urandom() < epsilon) {
        return urandom(0, n_actions);
    }
    return  MultinomialDistribution::generateInt(p);
}

int FixedContinuousPolicy::SelectAction(const Vector& next_state)
{
	state = next_state;	
	StatePolicy();
    if (epsilon_greedy && urandom() < epsilon) {
        return urandom(0, n_actions);
    }
	return  MultinomialDistribution::generateInt(p);
}

void FixedContinuousPolicy::Observe(const Vector& previous_state, const int& action, real r, const Vector& next_state)
{
    state = next_state;
	StatePolicy();
}

void FixedContinuousPolicy::Observe(real r, const Vector& next_state)
{
    state = next_state;
	StatePolicy();
}

void FixedContinuousPolicy::Reset(const Vector& start_state)
{
    state = start_state;
	StatePolicy();
}

real FixedContinuousPolicy::getActionProbability(const int& action) const
{
	assert(action >= 0 && action < n_actions);
	return p[action];
}

real FixedContinuousPolicy::getActionProbability(const Vector& start_state, const int& action) 
{
	assert(action >= 0 && action < n_actions);
	state = start_state;
	StatePolicy();
	return p[action];
}

void FixedContinuousPolicy::Show()
{
	for(int i = 0; i<n_actions; ++i){
		printf("%f ", getActionProbability(i));
	}
	printf("\n");
}

void FixedContinuousPolicy::StatePolicy()
{
	Vector Q(n_actions);
	bfs.Evaluate(state);
	Vector Phi = bfs.F();
	for(int i = 0; i< n_actions; ++i)
	{
		Q[i] = weights[(bfs.size() + 1)*i];
		for(int j = 0; j<bfs.size(); ++j)
		{
			Q[i] += Phi[j]*weights[i*(bfs.size() + 1) + j + 1];
		}
	}
	p.Clear();
	p[ArgMax(Q)] = 1.0;
}


// ---  Softmax --- //


///Create an uniform policy with a given number of basis functions.
SoftmaxContinuousPolicy::SoftmaxContinuousPolicy(int n_dimension_, int n_actions_, BasisSet<Vector, int>& bfs_)
	: ContinuousPolicy(n_dimension_, n_actions_, bfs_ )
{
	logmsg("Initialising softmax continuous policy");
	Vector w(n_actions*(bfs_.size() + 1));
	weights = w;
	p.Resize(n_actions);
}

SoftmaxContinuousPolicy::SoftmaxContinuousPolicy(int n_dimension_, int n_actions_, BasisSet<Vector, int>& bfs_, const Vector& weights_)
	: ContinuousPolicy(n_dimension_, n_actions_, bfs_, weights_)
{
	logmsg("Initialising softmax continuous policy");
	p.Resize(n_actions);
}

/// Destructor. Does nothing at the moment.
SoftmaxContinuousPolicy::~SoftmaxContinuousPolicy()
{
}

int SoftmaxContinuousPolicy::SelectAction()
{
	StatePolicy();
    return  MultinomialDistribution::generateInt(p);
}

int SoftmaxContinuousPolicy::SelectAction(const Vector& next_state)
{
	state = next_state;	
	StatePolicy();
	return  MultinomialDistribution::generateInt(p);
}

void SoftmaxContinuousPolicy::Observe(const Vector& previous_state, const int& action, real r, const Vector& next_state)
{
    state = next_state;
	StatePolicy();
}

void SoftmaxContinuousPolicy::Observe(real r, const Vector& next_state)
{
    state = next_state;
	StatePolicy();
}

void SoftmaxContinuousPolicy::Reset(const Vector& start_state)
{
    state = start_state;
	StatePolicy();
}

real SoftmaxContinuousPolicy::getActionProbability(const int& action) const
{
	assert(action >= 0 && action < n_actions);
	return p[action];
}

real SoftmaxContinuousPolicy::getActionProbability(const Vector& start_state, const int& action) 
{
	assert(action >= 0 && action < n_actions);
	state = start_state;
	StatePolicy();
	return p[action];
}

void SoftmaxContinuousPolicy::Show()
{
	for(int i = 0; i<n_actions; ++i){
		printf("%f ", getActionProbability(i));
	}
	printf("\n");
}

void SoftmaxContinuousPolicy::StatePolicy()
{
	Vector Q(n_actions);
	bfs.Evaluate(state); // state should be 2 dimensional
	Vector Phi = bfs.F();
	for(int i = 0; i< n_actions; ++i)
	{
		Q[i] = weights[(bfs.size() + 1)*i];
		for(int j = 0; j<bfs.size(); ++j)
		{
			Q[i] += Phi[j]*weights[i*(bfs.size() + 1) + j + 1];
		}
	}
	p.Clear();
	p[ArgMax(Q)] = 1.0;
}


// The gradient is phi(s,a) - sum_b phi(s,b) p(b|s)
const Vector SoftmaxContinuousPolicy::GradientUpdate(const Vector& s, const int& a, const real U)
{
	//bfs.Evaluate(state); // called by StatePolicy
	Vector tmp = state;
	state = s;
	StatePolicy();
	state = tmp;
	Vector phi = bfs.F();
	Vector delta(phi.Size()*n_actions);
	for (int b=0; b<n_actions; ++b) {
		Vector features(phi.Size()*n_actions);
		features.Insert(phi, b * phi.Size());
		real p = getActionProbability(b);
		if (b==a) {
			delta += features * (1.0-p);
		} else {
			delta -= features * p;
		}
	}
	return delta;
}
