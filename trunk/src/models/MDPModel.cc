// -*- Mode: c++ -*-
// copyright (c) 2005-2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: MDPModel.c,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MDPModel.h"
#include "Distribution.h"
#include "Random.h"

#include <iostream>

DiscreteMDP* MDPModel::CreateMDP()
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions, NULL, NULL);
    for (int i=0; i<n_states; ++i) {
        for (int a=0; a<n_actions; ++a) {
            SingularDistribution* ER = new SingularDistribution(getExpectedReward(i, a));
            mdp->addRewardDistribution(i, a, ER);
            real sum_p = 0.0;
            for (int j=0; j<n_states; ++j) {
                real p = getTransitionProbability(i, a, j);
                if (p) {
                    mdp->setTransitionProbability(i, a, j, p);
                    sum_p += p;
                    //printf("p(s'=%d|s=%d, a=%d)=%f\n", j, i, a, p);
                }
            }
            if (fabs(sum_p - 1.0) > 0.001) {
                printf ("!!!!!!! ");
                printf ("sum_s' p(s'|s=%d, a=%d) = %f\n", i, a, sum_p);
            }
        }
    }
    mdp->Check();
    return mdp;
}

GradientDescentMDPModel::GradientDescentMDPModel (int n_states, int n_actions, Distribution* initial_transitions, Distribution* initial_rewards) 
    : MDPModel (n_states, n_actions)
{
    this->initial_transitions = initial_transitions;
    this->initial_rewards = initial_rewards;

    N = n_states * n_actions;

    P = new real* [N];
    for (int i=0; i<N; i++) {
        P[i] = new real [n_states];
    }
    R = new real [N];

    Reset();
}

void GradientDescentMDPModel::Reset()
{
    for (int i=0; i<N; i++) {
        for (int j=0; j<n_states; j++) {
            P[i][j] = initial_transitions->generate();
        }
    }

    for (int i=0; i<N; i++) {
        R [i] = initial_rewards->generate();
    } 
}

GradientDescentMDPModel::~GradientDescentMDPModel()
{
    for (int i=0; i<N; i++) {
        delete [] P[i];
    }
    delete [] P;
    delete [] R;
}


void GradientDescentMDPModel::AddTransition(int s, int a, real r, int s2)
{
    int ID = getID (s, a);
    real* Ps=P[ID];
        
        
    SMART_ASSERT(s2>=0 && s2<n_states);
    // update transition probabilities
    real sum = 0.0f;
    for (int i=0; i<n_states; i++) {
        if (i==s2) {
            Ps[i] += alpha * (1.0f - Ps[i]);
        } else {
            Ps[i] -= alpha * Ps[i];
        }
        sum += Ps[i];
    }
        

    if (fabs(sum-1.0f)>0.01f) {         // should not be necessary to normalise
        std::cerr << "vector sum "
                  << sum << " violates constraint, renormalising.\n";
        real isum = 1.0f/sum;
        for (int i=0; i<n_states; i++) {
            Ps[i] *= isum;
        }
    }

    // update rewards (maybe replace with a gaussian model?)
    R[ID] += alpha * (r - R[ID]);
}





void GradientDescentMDPModel::ShowModel()
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
                  << "] = " << R[i] << std::endl; 
    }
}

real GradientDescentMDPModel::GenerateReward (int s, int a) const
{
    int ID = getID (s, a);
    return R[ID];
}

int GradientDescentMDPModel::GenerateTransition (int s, int a) const
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

real GradientDescentMDPModel::getTransitionProbability (int s, int a, int s2) const
{
    int ID = getID (s, a);                
    assert (s2>=0 && s2<n_states);
    return P[ID][s2];
}

real GradientDescentMDPModel::getExpectedReward (int s, int a) const
{
    int ID = getID (s, a);
    return R[ID];
}
