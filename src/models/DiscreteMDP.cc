// -*- Mode: c++ -*-
// copyright (c) 2005-2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
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
      N(n_states * n_actions),
	  reward_distribution(n_states, n_actions),
      transition_distribution(n_states, n_states)
{   
	if (initial_transitions) {
		Serror("Not implemented\n");
	}
    //real p = 1.0 / (real) n_states;
    reward = 0.0;
    state = 0;
}

/** Partially copies an MDP.
    
    Since there is no way to clone the distribution pointers,
    we actuallcd al create new, singular distributions instead.
*/
DiscreteMDP::MDP(const MDP<int,int>& mdp)
    : n_states(mdp.n_states),
      n_actions(mdp.n_actions),
	  reward_distribution(n_states, n_actions),
	  transition_distribution(mdp.transition_distribution)
{
    reward = 0.0;
	state = 0;
	reward_distribution = mdp.reward_distribution;
}


DiscreteMDP::MDP (const std::vector<const MDP<int,int>*> &mdp_list,
                  const Vector& w)
    : n_states(mdp_list[0]->n_states),
      n_actions(mdp_list[0]->n_actions),
      reward_distribution(n_states, n_actions),
	  transition_distribution(n_states, n_actions)
{
    int n_mdps = mdp_list.size();
    for (int i=0; i<n_mdps; ++i) {
        if (n_states != mdp_list[i]->n_states) {
            Serror("Number of states don't match\n");
            exit(-1);
        }
        if (n_actions != mdp_list[i]->n_actions) {
            Serror("Number of actions don't match\n");
            exit(-1);
        }
    }

    assert(w.Sum() == 1.0);

    for (int i=0; i<n_states; ++i) {
        for (int a=0; a<n_actions; ++a) {
            for (int j=0; j<n_states; ++j) {
                real p_ija = 0.0;
                for (int m=0; m<n_mdps; ++m) {
                    p_ija += w(m) * mdp_list[m]->getTransitionProbability(i, a, j);
                }
                setTransitionProbability(i, a, j, p_ija);
            }
            real r_ia = 0.0;
            for (int m=0; m<n_mdps; ++m) {
                r_ia += w(m) * mdp_list[m]->getExpectedReward(i, a);
            }
            setFixedReward(i, a, r_ia);
        }
    }    

    reward =0.0;
    state = 0;
}

DiscreteMDP::~MDP()
{
}


/** From Putterman, 8.5.4.
 
    The gain value function is now proportional to tau.
*/
void DiscreteMDP::AperiodicityTransform(real tau)
{
    //DISABLED_ASSERT(tau > 0 && tau < 1)(tau);
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
			reward_distribution.setFixedReward(s, a, tau * reward_distribution.expected(s, a));
		}
	}
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            for (int j=0; j<n_states; j++) {
                real delta = 0;
                if (j==s) {
                    delta = 1.0;
                }
				real P = (1-tau)*delta + tau*getTransitionProbability(s,a,j);
                transition_distribution.SetTransition(s,a,j,P);
													  
            }
        }
    }
}


void DiscreteMDP::ShowModel() const
{
	real threshold = 10e-5;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            std::cout << "(" << s << "," << a << ") -> ";
			real sum = 0.0;
            for (int j=0; j<n_states; j++) {
                real p = getTransitionProbability(s,a,j);
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
            //int i = getID(s,a);
            for (int j=0; j<n_states; j++) {
                real p = getTransitionProbability(s,a,j);
                if (p>0.000001) {
                    fprintf (fout,
                             "s%d -> s%d [label = \" p=%.2f, r=%.1f\", color=%s];\n",
                             s, j, p, reward_distribution.expected(s,a), colour[colour_id]);
                }
            }
        }
    }

}



real DiscreteMDP::generateReward (const int& s, const int& a) const
{
    return reward_distribution.generate(s, a);
}

int DiscreteMDP::generateState (const int& s, const int& a) const
{
    //int ID = getID (s,a);
    real sum = 0.0f;
    real X = urandom();

    int select = 0;
    for (int i=0; i<n_states; i++) {
        sum += getTransitionProbability(s, a, i);
        if (X<sum) {
            select = i;
            break;
        }
    }
    return select;
}


#if 0
bool DiscreteMDP::Check() const
{
    return true;
}
#else
bool DiscreteMDP::Check() const
{
    real threshold = 0.001;
    bool flag = true;
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            real sum = 0.0;
            for (int s2=0; s2<n_states; s2++) {
                real p = getTransitionProbability(s, a, s2);
                //printf ("%d %d : - %f -> %d\n", s, a, p, s2);
                sum += p;
            }
            assert(fabs(sum - 1.0f) <= threshold);//(s)(a)(sum);
            if (fabs(sum - 1.0f) > threshold) {
                flag = false;
				Serror("transition s:%d a:%d = %f\n", s, a, sum);
            }
        }
    }
    logmsg("MDP Check, %d states, %d actions: %d\n", n_states, n_actions, flag);
    return flag;
}
#endif
