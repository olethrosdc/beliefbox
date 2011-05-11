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
#include "Matrix.h"

#include <iostream>


DiscreteMDP::MDP (int n_states_, int n_actions_, real** initial_transitions)
	: n_states(n_states_),
      n_actions(n_actions_),
      P(n_states*n_actions, n_states),
      next_states(n_states * n_actions),
      N(n_states * n_actions),
	  reward_distribution(n_states, n_actions)
{   
    next_states.resize(N);
    
    //real p = 1.0 / (real) n_states;
    if (initial_transitions) {
        for (int i=0; i<N; i++) {
            //P[i] = &P_data[i*n_states];
            for (int j=0; j<n_states; j++) {
                //P[i][j] = initial_transitions[i][j];
				P(i,j) = initial_transitions[i][j];
            }
        } 
    } else {
        int i=0;
        for (int s=0; s<n_states; s++) {
            for (int a=0; a<n_actions; a++, i++) {
                for (int j=0; j<n_states; j++) {
						P(i,j) = 0.0;
                }
            }
        }
    }

#if 0    
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
#endif
}

/** Partially copies an MDP.
    
    Since there is no way to clone the distribution pointers,
    we actually create new, singular distributions instead.
*/
DiscreteMDP::MDP(const MDP<int,int>& mdp)
    : n_states(mdp.n_states),
      n_actions(mdp.n_actions),
	  reward_distribution(n_states, n_actions)
{
    N = n_states * n_actions;

	state = 0;
    next_states.resize(N);
	P = mdp.P;
	reward_distribution = mdp.reward_distribution;
}


DiscreteMDP::~MDP()
{
}


/** From Putterman, 8.5.4.
 
    The gain value function is now proportional to tau.
 */
void DiscreteMDP::AperiodicityTransform(real tau)
{
    SMART_ASSERT(tau > 0 && tau < 1)(tau);
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
			reward_distribution.setFixedReward(s, a, tau * reward_distribution.expected(s, a));
		}
	}
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int i = getID(s,a);
            for (int j=0; j<n_states; j++) {
                real delta = 0;
                if (j==s) {
                    delta = 1.0;
                }
                P(i,j) = (1-tau)*delta + tau*P(i,j);
            }
        }
    }
}

void DiscreteMDP::ShowModel() const
{
	real threshold = 0.001;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            int i = getID(s,a);
            std::cout << "(" << s << "," << a << ") -> ";
			real sum = 0.0;
            for (int j=0; j<n_states; j++) {
                real p = P(i,j);
				sum += p;
                if (p>threshold) {
                    std::cout << j << " (" << p << ") ";
                }
            }
			if (fabs(sum - 1.0) > threshold) {
				std::cout << "# ERR";
			}
            std::cout << std::endl;
        }
    }
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            //int i = getID(s,a);
            std::cout << "R[" << s << "," << a << "] = "
                      << reward_distribution.expected(s,a) << std::endl; 
        }
    }
}

void DiscreteMDP::dotModel(FILE* fout) const
{
    const char* colour[] ={
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
                real p = P(i, j);
                if (p>0.000001) {
                    fprintf (fout,
                             "s%d -> s%d [label = \" p=%.2f, r=%.1f\", color=%s];\n",
                             s, j, p, reward_distribution.expected(s,a), colour[colour_id]);
                }
            }
        }
    }

}



real DiscreteMDP::generateReward (int s, int a) const
{
    return reward_distribution.generate(s, a);
}

int DiscreteMDP::generateState (int s, int a) const
{
    int ID = getID (s,a);
    real sum = 0.0f;
    real X = urandom();

    int select = 0;
    for (int i=0; i<n_states; i++) {
        sum += P(ID, i);
        if (X<sum) {
            select = i;
            break;
        }
    }
    return select;
}


#ifdef NDEBUG
bool DiscreteMDP::Check() const
{
    return true;
}
#else
bool DiscreteMDP::Check() const
{
    real threshold = 0.001;
    bool flag = false;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            real sum = 0.0;
            for (int s2=0; s2<n_states; s2++) {
                real p = getTransitionProbability(s, a, s2);
                sum += p;
            }
            assert(fabs(sum - 1.0f) <= threshold);//(s)(a)(sum);
            if (fabs(sum - 1.0f) > threshold) {
                flag = true;
            }
        }
    }
    return flag;
}
#endif
