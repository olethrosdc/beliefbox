// -*- Mode: c++ -*-
// copyright (c) 2005-2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/************************************************************************
 *                                                                      *
 * This program is free software; you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation; either version 2 of the License, or    *
 * (at your option) any later version.                                  *
 *                                                                      *
 ************************************************************************/

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
    
    P.resize(N);//= new real* [N];
    P_data.resize(N*n_states);
    state = 0;
    
    next_states.resize(N);
    
    for (int i=0; i<N; i++) {
        P[i] = &P_data[i*n_states];
        if (initial_transitions) {
            for (int j=0; j<n_states; j++) {
                P[i][j] = initial_transitions[i][j];
            }
        } else {
            for (int j=0; j<n_states; j++) {
                P[i][j] = 0.0;//1.0 / (real) n_actions;
            }
        }
    }
    
    ER.resize(N);
    R.resize(N);
    if (initial_rewards) {
        for (int i=0; i<N; i++) {
            R[i] = initial_rewards[i];
            ER[i] = R[i]->getMean();
        }
    } else {
        for (int i=0; i<N; i++) {
            R[i] = NULL;
            ER[i] = 0xBADFEED;
        }
    }
}


DiscreteMDP::~MDP()
{
}



void DiscreteMDP::ShowModel() const
{
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int i = getID(s,a);
            std::cout << "(" << s << "," << a << ") :";
            for (int j=0; j<n_states; j++) {
                real p = P[i][j];
                if (p>0.01) {
                    std::cout << j << " (" << p << ") ";
                }
            }
            std::cout << std::endl;
        }
    }
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int i = getID(s,a);
            std::cout << "R[" << s << "," << a << "] = "
                      << ER[i] << std::endl; 
        }
    }
}

void DiscreteMDP::dotModel(FILE* fout) const
{
    char* colour[] ={
        "red",
        "green",
        "blue",
        "yellow",
        "magenta",
        "cyan",
        "black"};

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int colour_id = a % 7;
            int i = getID(s,a);
            for (int j=0; j<n_states; j++) {
                real p = P[i][j];
                if (p>0.000001) {
                    fprintf (fout,
                             "s%d -> s%d [label = \" p=%.2f, r=%.1f\", color=%s];\n",
                             s, j, p, ER[i], colour[colour_id]);
                }
            }
        }
    }

}



real DiscreteMDP::generateReward (int s, int a) const
{
    int ID = getID (s, a);
    assert (R[ID]);
    //dprint ("ID: %d : %p\n", ID, R[ID]);
    return ER[ID];
    //    return R[ID]->generate();
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


#ifdef NDEBUG
void DiscreteMDP::Check() const
{
    
}
#else
void DiscreteMDP::Check() const
{
    real threshold = 0.001;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            real sum = 0.0;
            for (int s2=0; s2<n_states; s2++) {
                real p = getTransitionProbability(s, a, s2);
                sum += p;
            }
            SMART_ASSERT(fabs(sum - 1.0f) <= threshold)(s)(a)(sum);
        }
    }
}
#endif
