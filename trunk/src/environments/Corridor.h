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

#ifndef CORRIDOR_H
#define CORRIDOR_H

#include "real.h"

class DiscretePOMDP;
class RandomNumberGenerator;

class Corridor
{
protected:
    int n_states;
    int observation;  ///< current observation
    int state; ///< current state
    real reward; ///< current reward
    real randomness;
    RandomNumberGenerator* rng;
    DiscretePOMDP* pomdp;
public:
    Corridor(int n_states_, real randomness_, RandomNumberGenerator* rng_)
        : n_states(n_states_),
          randomness(randomness_),
          rng(rng_)
    {
        printf ("# Making Corridor of length %d\n", n_states);
        assert(n_states > 0);
        pomdp = new DiscretePOMDP(n_states, 2, 2);
        int n_obs = 2;
        // clear all
        for (int s=0; s<n_states; ++s) {
            for (int s2=0; s2<n_states; ++s2)  {
                pomdp->setNextStateProbability(s, 0, s2, 0);
                pomdp->setNextStateProbability(s, 1, s2, 0);
            }
            for (int x=0; x<n_obs; ++x)  {
                pomdp->setObservationProbability(s, 0, x, 0);
                pomdp->setObservationProbability(s, 0, x, 0);
            }
        }
        
        // set states
        for (int s=0; s<n_states-1; ++s) {
            pomdp->setNextStateProbability(s+1, 0, s, 1);
            pomdp->setNextStateProbability(s, 1, s+1, 1);
        }
        pomdp->setNextStateProbability(0, 0, 0, 1);
        pomdp->setNextStateProbability(n_states-1, 1, n_states-1, 1);

        // set obs
        for (int s=0; s<n_states-1; ++s) {
            // see nothing when no obstacle
            pomdp->setObservationProbability(s+1, 0, 0, 1 - randomness);
            pomdp->setObservationProbability(s,   1, 0, 1 - randomness);
            pomdp->setObservationProbability(s+1, 0, 1, randomness);
            pomdp->setObservationProbability(s,   1, 1, randomness);
        }
        // see something when obstacle
        pomdp->setObservationProbability(0, 0, 1, 1 - randomness);
        pomdp->setObservationProbability(0, 0, 0, randomness);
        pomdp->setObservationProbability(n_states-1, 1, 1, 1 - randomness);
        pomdp->setObservationProbability(n_states-1,  1, 0, randomness);

        pomdp->check();
    }
    ~Corridor()
    {
        delete pomdp;
    }
    void Reset()
    {
        state = 0;
        observation = 0;
    }
    DiscretePOMDP* getPOMDP()
    {
        return pomdp;
    }
    int getObservation()
    {
        return observation;
    }
    real getReward()
    {
        return reward;
    }
    int getNStates()
    {
        return n_states;
    }
    int getNActions()
    {
        return 2;
    }
    int getNObservations()
    {
        return 2;
    }
    // give a reward of 1 when you are at the last state
    bool Act(int action)
    {
        switch(action) {
        case 0:
            if (rng->uniform() < randomness) {
                state++;
            } else {
                state--;
            } 
            break;
        case 1:
            if (rng->uniform() < randomness) {
                state--;
            } else {
                state++;
            }
            break;
        default:
            break;
        }

        observation = 0;
        reward = rng->uniform()<0.5;
        if (state <= 0) {
            state = 0;
            observation = 1;
            reward = 0;
        } else if (state >= n_states - 1) {
            state = n_states - 1;
            observation = 1;
        }
        
        if (state == n_states - 1) {
            reward = 1;
        }
            
        if (rng->uniform() < randomness) {
            observation = 1 - observation;
        }
        pomdp->setObservation(observation);
        pomdp->setReward(reward);
        pomdp->setState(state);
        return 0;
    }
};

#endif
