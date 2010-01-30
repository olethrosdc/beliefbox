/* -*- Mode: c++;  -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef POMDP_BELIEF_PREDICTOR_H
#define POMDP_BELIEF_PREDICTOR_H

#include "POMDPBeliefState.h"
#include "DiscretePOMDP.h"
#include "FactoredPredictor.h"

/// Abstract class for prediction with actios
class POMDPBeliefPredictor : public FactoredPredictor
{
    
protected:
    DiscretePOMDP* pomdp;
    DiscretePOMDPBeliefState* belief_state;
public:
    POMDPBeliefPredictor(DiscretePOMDP* pomdp_)
        : pomdp(pomdp_)
    {
        belief_state = new DiscretePOMDPBeliefState(pomdp);
    }
    virtual ~POMDPBeliefPredictor()
    {
        delete belief_state;
    }
    
    /* Training and generation */
    virtual real Observe (int prd)
    {
        return 1.0 / (real) pomdp->getNObservations();
    }
    virtual real Observe (int act, int prd)
    {
        return belief_state->Observe(act, prd, 0);
    }
    virtual real ObservationProbability (int act, int x) 
    {
        return belief_state->ObservationProbability(act, x, 0);
    }
    virtual void Reset()
    {                                          
        belief_state->Reset();
    }
    
}; 



#endif
