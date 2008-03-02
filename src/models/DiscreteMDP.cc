// -*- Mode: c++ -*-
// copyright (c) 2005-2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMDP.h"
#include "Distribution.h"
#include "Random.h"
#include "SmartAssert.h"
#include <iostream>

DiscreteMDP::MDP (int n_states, int n_actions, real** initial_transitions, Distribution** initial_rewards) 
{   
    this->n_states = n_states;
    this->n_actions = n_actions;
    
    N = n_states * n_actions;

    P = new real* [N];
    P_data = new real[N*n_states];
	state = 0;

    for (int i=0; i<N; i++) {
        P[i] = &P_data[i*n_states];
        if (initial_transitions) {
            for (int j=0; j<n_states; j++) {
                P[i][j] = initial_transitions[i][j];
            }
        } else {
            for (int j=0; j<n_states; j++) {
                P[i][j] = 1.0 / (real) n_actions;
            }
        }
    }

    ER = new real [N];
    R = new Distribution* [N];
    if (initial_rewards) {
        for (int i=0; i<N; i++) {
            R[i] = initial_rewards[i];
            ER[i] = R[i]->getMean();
            //printf ("%d %f\n", i, ER[i]);
        }
    }

}


DiscreteMDP::~MDP()
{
    delete [] P_data;
    delete [] P;
    delete [] R;
	delete [] ER;
}


void DiscreteMDP::setTransitionProbability(int s, int a, int s2, real p)
{
    int ID = getID (s, a);
    real* Ps=P[ID];
    SMART_ASSERT(s2>=0 && s2<n_states)(s2);
    Ps[s2] = p;
}        

void DiscreteMDP::setRewardDistribution(int s, int a, Distribution* reward)
{   
    int ID = getID (s, a);
    R[ID] = reward;
    ER[ID] = reward->getMean();
    //printf ("%d %f\n", ID, ER[ID]);
}

void DiscreteMDP::ShowModel() const
{
    for (int i=0; i<N; i++) {
        std::cout << i << ":";
        for (int j=0; j<n_states; j++) {
            real p = P[i][j];
            if (p<0.01) p =0.0f;
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }

    for (int i=0; i<N; i++) {
        std::cout << "R[" << i
                  << "] = " << ER[i] << std::endl; 
    }
}

real DiscreteMDP::generateReward (int s, int a) const
{
    int ID = getID (s, a);
    return R[ID]->generate();
}

int DiscreteMDP::generateState (int s, int a) const
{
    int ID = getID (s,a);
    real* Ps=P[ID];
    real sum = 0.0f;
    real X = urandom();

    int select = 0;
    for (int i=0; i<n_states; i++) {
        sum += Ps[i];
        if (X<sum) {
            select = i;
            break;
        }
    }
    return select;
}

real DiscreteMDP::getTransitionProbability (int s, int a, int s2) const
{
    int ID = getID (s, a);                
    assert (s2>=0 && s2<n_states);
    return P[ID][s2];
}

real DiscreteMDP::getExpectedReward (int s, int a) const
{
    int ID = getID (s, a);
    return ER[ID];
}

void DiscreteMDP::Check() const
{
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            real sum = 0.0;
			printf ("E[r|s=%d, a=%d] = %f\n", s, a, getExpectedReward(s, a));
            for (int s2=0; s2<n_states; s2++) {
                real p = getTransitionProbability(s, a, s2);
                printf ("P[s'=%d| s=%d, a=%d] = %f\n", s2, s, a, p);
                sum += p;
            }
            SMART_ASSERT(fabs(sum - 1.0f) <= 0.001)(s)(a)(sum);
        }
    }
}
