// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef POMDP_H
#define POMDP_H

#include "real.h"
#include "SmartAssert.h"
#include "MDP.h"


template<typename ObservationType, typename StateType, typename ActionType>
class POMDP : public MDP<StateType, ActionType> {
protected:
	ObservationDistribution<StateType, ObservationType>& observation_distribution;
public:
    POMDP(TransitionDistribution<StateType, ActionType>& transition_distribution_,
		  RewardDistribution<StateType, ActionType>& reward_distribution_,
		  ObservationDistribution<StateType, ObservationType>& observation_distribution_)
        : MDP<StateType, ActionType> (transition_distribution_, reward_distribution),
		  observation_distribution(observation_distribution_)

		  virtual ~POMDP()
    {
    }
};

