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

#ifndef POMDP_WRAPPER_H
#define POMDP_WRAPPER_H

#include "Environment.h"

/** Wrap a POMDP to look like an MDP to a normal agent.
    
    Just map observations to states.
 */
template <typename S, typename A, typename P>
class POMDPWrapper : public Environment<S, A>
{
public:
    P& pomdp;
    POMDPWrapper(P& pomdp_) 
        : pomdp(pomdp_)
    {
        n_states = pomdp.getNObs();
        n_actions = pomdp.getNActions();
        state = pomdp.getObservation();
        reward = pomdp.getReward();
    }
    virtual ~POMDPWrapper()
    {
    }
    /// put the environment in its "natural: state
    virtual void Reset() 
    {
        pomdp.Reset();
        state = pomdp.getObservation();
        reward = pomdp.getReward();
        S x = getState();
    }

    /// returns true if the action succeeds, false if it does not.
    ///
    /// The usual of false is that the environment is in a terminal
    /// absorbing state.
    virtual bool Act(A action) 
    {
        bool act = POMDP.Act(action);
        state = pomdp.getObservation();
        reward = pomdp.getReward();
    }
    

    virtual const char* Name()
    {
        return "POMDP Wrapper";
    }

};

#endif
