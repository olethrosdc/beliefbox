// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: MDP.h,v 1.3 2006/11/06 23:42:32 olethros Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TRANSITION_DISTRIBUTION_H
#define TRANSITION_DISTRIBUTION_H

#include "DiscreteStateSet.h"

template <typename StateType, typename ActionType>
class TransitionDistribution
{
 public:
    virtual ~TransitionDistribution() {}
    virtual StateType generate(StateType state, ActionType action) = 0;
    virtual real pdf(StateType state, ActionType action, StateType next_state) = 0;
};

#if 0
class TransitionDistribution<int, int>
{
public:
  int n_states;
  int n_actions;
  real epsilon;
  DiscreteStateSet states;
  TransitionDistribution(int n_states_, int n_actions_, real epsilon_)
    : n_states(n_states_),
      n_actions(n_actions_),
      epsilon(epsilon)
  {
  }
  virtual ~TransitionDistribution() {}
  virtual StateType generate(int state, int action) 
  {
  }
  virtual real pdf(int state, int action, int next_state)
  {
    DiscreteStateSetRef
  }
};
#endif
#endif
