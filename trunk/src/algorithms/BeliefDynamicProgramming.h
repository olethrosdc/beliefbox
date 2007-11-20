// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BELIEF_DYNAMIC_PROGRAMMING_H
#define BELIEF_DYNAMIC_PROGRAMMING_H

/** Dynamic programming using beliefs.

We start with an initial belief \f$b_0\f$ which is a distribution over
possible worlds.  Then, we expand the tree according to the possible
outcomes at each time \f$t\f$ to obtain \f$\{b_i(t)\}\f$.

The question is, how do we expand the tree?  We have to model all the
variables as observations.  If we are talking about reinforcement
learning, then those variables will be the tuple \f$(S_t, A_t,
S_{t+1}, R_t)\f$.  It is necessary that a distribution be given for
these tuples.

*/
template <typename StateType, typename ActionType>
class BeliefDynamicProgramming
{
protected:
	Belief* initial_belief;
public:
	BeliefDynamicProgramming(Belief* initial_belief);
	~virtual BeliefDynamicProgramming(Belief* initial_belief);
};

/** Belief MDPs for Continuous Bandit Problems with Bernoulli rewards 
	Assumes that the space is n-dimensional */
class BDP_ContinuousBandit
{
protected:
	int n_dimensions; ///< number of dimensions
public:
	BDP_ContinuousBandit(int d) 
	{
		n_dimensions = d;
	}
	virtual ~BDP_ContinuousBandit()
	{}
	virtual void Observe(Vector a, real r);
	
};


#endif

