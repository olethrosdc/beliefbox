// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REWARD_BELIEF_H
#define REWARD_BELIEF_H

/** A belief about rewards in discrete spaces.
    
    We wish to split this belief into bits.
 */
class DiscreteRewardBelief
{
protected:
	int n_states; ///< number of states
	int n_actions; ///< number of actions
    std::vector<Distribution*> R; ///< reward distribution
	std::vector<Distribution*> distribution_vector; ///< for malloc
	Vector ER; ///< expected reward
    inline int getID (int s, int a) const
    {
        
        assert(s>=0 && s<n_states);
        assert(a>=0 && a<n_actions);
        return s*n_actions + a;
    }

public:
    // only set up
    DiscreteRewardBelief(int n_states_, int n_actions_)
        : n_states(n_states_), n_actions(n_actions_)
    {}
        
    virtual ~DiscreteRewardBelief();
    virtual real generate(int state, int action) const;
    virtual real expected(int state, int action) const;
    virtual real pdf(int state, int action, real reward) const;
	void setRewardDistribution(int s, int a, Distribution* reward);
	void addRewardDistribution(int s, int a, Distribution* reward);
	void addFixedReward(int s, int a, real reward);
	void setFixedReward(int s, int a, real reward);
};
