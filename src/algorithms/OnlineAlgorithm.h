// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/**
   \defgroup ReinforcementLearning Reinforcement learning algorithms.

   This group contains various algorithms for the reinforcement
   learning problem. These algorithms usually employ models,
   optimisation methods and heuristics in order to learn how to act
   near-optimally.

*/
/* @{ */

#ifndef ONLINE_ALGORITHM_H
#define ONLINE_ALGORITHM_H

/** Online algorithm template.

    This is simply a template for reinforcement learning algorithms.
    This acts as a base class for all algorithms that are going to be
    solving reinforcement learning problems. That means that given a
    sequence of observations and rewards, the produce a sequence of
    actions. Other algorithms related to reinforcement learning, such
    as value iteration, are not derived from this class, but can be
    used in the implementation of any derived class.  */
template <typename A, typename S>
class OnlineAlgorithm
{
public:
    OnlineAlgorithm()
    {
    }
    virtual ~OnlineAlgorithm()
    {
    }
    /// call this at the end of an episode.
    virtual void Reset()
    {
    }
    /// Partial observation (when another algorithm is used to take the action, or as a helper function)
    virtual real Observe (real reward, const S& next_state, const A& next_action) = 0;
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual A Act(real reward, const S& next_state) = 0;
	/// Some algorithms may implement getting the value of a state and action. Implementation is optional.
    virtual real getValue (const S& state, const A& action)
	{
		Swarning("not implemented\n");
		return 0;
	}
    /// Some algorithms may implement a different strategy when the reward matrix SxA is given. Implementation is not obligatory.
    virtual void setFixedRewards(const Matrix& rewards) 
    {
        Swarning("not implemented\n");
    }
};

/// @}

#endif
