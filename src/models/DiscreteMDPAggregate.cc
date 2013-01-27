// -*- Mode: c++ -*-
// copyright (c) 2005-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMDPAggregate.h"
#include "Random.h"
#include "Gridworld.h"


/// This constructor makes a random set 
DiscreteMDPAggregate::DiscreteMDPAggregate (int n_aggregated_states, int n_states, int n_actions, int init_transition_count) :
    DiscreteMDPCounts(n_states, n_actions, init_transition_count)

{
    mdp_dbg("Creating DiscreteMDPAggregate from %d states to %d sets and %d actions\n",  n_aggregated_states, n_states, n_actions);
    this->n_aggregated_states = n_aggregated_states;
    X.resize(n_states);
    state_map.resize(n_aggregated_states);

    BuildRandomAggregate();
}

DiscreteMDPAggregate::DiscreteMDPAggregate (Gridworld& gridworld, int n_aggregated_states, int n_states, int n_actions, int init_transition_count) 
 : DiscreteMDPCounts(n_states, n_actions, init_transition_count)
{
    mdp_dbg("Creating DiscreteMDPAggregate from %d states to %d sets and %d actions\n",  n_aggregated_states, n_states, n_actions);
    this->n_aggregated_states = n_aggregated_states;
    X.resize(n_states);
    state_map.resize(n_aggregated_states);

    BuildGridworldAggregate(gridworld);
}

void DiscreteMDPAggregate::BuildRandomAggregate()
{
    for (int i=0; i<n_aggregated_states; i++) {
        int zeros = 0;
        for (int x=0; x<n_states; x++) {
            if (X[x].size()==0) {
                zeros++;
            }
        }
        int s = (int) floor(urandom() * ((real) n_states));
        while (zeros && X[s].size()!=0) {
            s = (s+1) % n_states;
        }
        //DISABLED_ASSERT(s>=0 && s<n_states);
            
        X[s].add(i);
        state_map[i] = s;
    }

}

void DiscreteMDPAggregate::BuildGridworldAggregate(Gridworld& gridworld)
{
    int width = gridworld.getWidth();
    int height = gridworld.getHeight();

    // keep one state for the termination and other states.
    int n_divisions_x = (int) floor(sqrt(n_states-1));//width*height));
    int n_divisions_y = (n_states-1)/n_divisions_x;
    
    int remaining_states = n_states - n_divisions_x*n_divisions_y;
    printf("divx: %d, divy: %d, remaining: %d/%d\n",
           n_divisions_x,
           n_divisions_y,
           remaining_states,
           n_states);
    for (int i=0; i<n_aggregated_states; i++) {
        state_map[i] = -1;
    }
    int x_dist = width / n_divisions_x;
    int y_dist = height / n_divisions_y;
    for (int y=0; y<n_divisions_y; ++y) {
        for (int x=0; x<n_divisions_x; ++x) {
            for (int y_orig=y*y_dist;
                 y_orig < std::min(height,(y+1)*y_dist);
                 ++y_orig) {
                for (int x_orig=x*x_dist;
                     x_orig < std::min(width,(x+1)*x_dist);
                     ++x_orig) {
                    int s = y*n_divisions_x + x;
                    int i = gridworld.getState(x_orig, y_orig);
                    X[s].add(i);
                    state_map[i] = s;
                }
            }
        }
    }
    for (int i=0; i<n_aggregated_states; i++) {
        // only fill in unmapped states
        if (state_map[i] >= 0) {
            continue;
        }
        int zeros = 0;
        for (int x=0; x<n_states; x++) {
            if (X[x].size()==0) {
                zeros++;
            }
        }
        int s = (int) floor(urandom() * ((real) n_states));
        // make sure you fill in at least one empty aggregated state
        while (zeros && X[s].size()!=0) {
            s = (s+1) % n_states;
        }
        //DISABLED_ASSERT(s>=0 && s<n_states);
        X[s].add(i);
        state_map[i] = s;
    }
    
}

DiscreteMDPAggregate::~DiscreteMDPAggregate()
{
}

void DiscreteMDPAggregate::AddTransition(int s, int a, real r, int s2)
{
    DiscreteMDPCounts::AddTransition(Aggregate(s), a, r, Aggregate(s2));
}

//void DiscreteMDPAggregate::SetNextReward(int s, int a, real r)
//{
//    DiscreteMDPCounts::SetNextReward(Aggregate(s), a, r);
//}

real DiscreteMDPAggregate::GenerateReward (int s, int a) const
{
    return DiscreteMDPCounts::GenerateReward(Aggregate(s), a);
}

int DiscreteMDPAggregate::GenerateTransition (int s, int a) const
{
    return DiscreteMDPCounts::GenerateTransition(Aggregate(s), a);
}

real DiscreteMDPAggregate::getTransitionProbability (int s, int a, int s2) const
{
    int x = Aggregate(s);
    int x2 = Aggregate(s2);
    int n = X[x2].size();
    real p = 0.0;
    if (n > 0) {
        p = DiscreteMDPCounts::getTransitionProbability(x, a, x2) / (real) n; 
    }
    mdp_dbg ("\n P(s'=%d | s=%d, a=%d) = (1/%d) P(x'=%d | x=%d, a=%d) = %f\n",
             s2, s, a,
             n, x2, x, a,
             p);
    return p;
}

Vector DiscreteMDPAggregate::getTransitionProbabilities (int s, int a) const
{
    return DiscreteMDPCounts::getTransitionProbabilities(Aggregate(s), a);
}

real DiscreteMDPAggregate::getExpectedReward (int s, int a) const
{
    return DiscreteMDPCounts::getExpectedReward(Aggregate(s), a);
}

void DiscreteMDPAggregate::Reset()
{
    DiscreteMDPCounts::Reset();
}
