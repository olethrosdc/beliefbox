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
        : Environment<S, A>(pomdp_.getNObs(), pomdp_.getNActions()),
          pomdp(pomdp_)
    {
        Environment<S,A>::state = pomdp.getObservation();
        Environment<S,A>::reward = pomdp.getReward();
    }
    virtual ~POMDPWrapper()
    {
    }
    /// put the environment in its "natural: state
    virtual void Reset() 
    {
        pomdp.Reset();
        Environment<S,A>::state = pomdp.getObservation();
        Environment<S,A>::reward = pomdp.getReward();
        S x = Environment<S,A>::getState();
    }

    /// returns true if the action succeeds, false if it does not.
    ///
    /// The usual of false is that the environment is in a terminal
    /// absorbing state.
    virtual bool Act(A action) 
    {
        bool act = pomdp.Act(action);
        Environment<S,A>::state = pomdp.getObservation();
        Environment<S,A>::reward = pomdp.getReward();
        return act;
    }
    

    virtual const char* Name()
    {
        return "POMDP Wrapper";
    }

};

#endif
