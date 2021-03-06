/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "PolicyEvaluation.h"
#include "ValueIteration.h"
#include "AverageValueIteration.h"
#include "Gridworld.h"
#include "InventoryManagement.h"
#include "DiscreteChain.h"
#include "RandomMDP.h"
#include "DiscretePolicy.h"
#include "MersenneTwister.h"
#include "MDPWrapper.h"

int main (void)
{
    int period = 20;
    int max_items = 5;
    real gamma = 0.95;
    real demand = 0.1;
    real random = 0.0;
    real pit = -100.0;
    real goal = 0.0;
    real step = -0.1;
	real margin = 1.1;
    real epsilon = 1e-6;

    MersenneTwisterRNG rng;
    InventoryManagement inventory_management (period, max_items, demand, margin);

    Gridworld grid_world("/home/olethros/projects/beliefbox/dat/maze2", random, pit, goal, step);
    DiscreteChain chain(5);
    RandomMDP random_mdp(2, 8, 0.0, 0.1, 0, 1, &rng, false);
    
    //DiscreteMDP* mdp = random_mdp.getMDP();
    //DiscreteMDP* mdp = inventory_management.getMDP();
    DiscreteMDP* mdp = grid_world.getMDP();
    //DiscreteMDP* mdp = chain.getMDP();

    int n_states = mdp->getNStates();
    int n_actions = mdp->getNActions();
    
    std::vector<Vector> p(n_states);
    for (int s=0; s<n_states; s++) {
        p[s].Resize(n_actions);
        for (int a=0; a<n_actions; a++) {
            p[s][a] = 1.0 / (real) n_actions;
	    //printf ("R(%d,%d) = %f\n", s, a, mdp->getExpectedReward(s, a));
        }
    }

    FixedDiscretePolicy* policy = new FixedDiscretePolicy(p);

    PolicyEvaluation policy_evaluation(policy, mdp, gamma);

    ValueIteration value_iteration(mdp, gamma);
    //AverageValueIteration value_iteration(mdp, false, false);

    std::vector<real> Q(n_actions);

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q[a] = policy_evaluation.getValue(s,a);
        }
        int a_max = ArgMax(Q);
        std::cout << s << " : " << policy_evaluation.getValue(s)
                  << " | "
                  << a_max << ":" << Q[a_max]
                  << std::endl;
    }

#if 1
    
   for (int iter=0; iter < 1; iter++) {
       printf ("ITER: %d, V = ", iter);
       value_iteration.ComputeStateValues(epsilon, -1);
       FixedDiscretePolicy* vi_policy = value_iteration.getPolicy();
       policy_evaluation.SetPolicy(vi_policy);
       policy_evaluation.ComputeStateValues(epsilon, -1);
       for (int s=0; s<mdp->getNStates(); s++) {
           printf ("%.1f (%.1f) ",
                   value_iteration.getValue(s),
                   policy_evaluation.getValue(s));
       }
       delete vi_policy;

       printf(", D = %.1f, b = %.1f\n",
              value_iteration.Delta,
              value_iteration.baseline);
   }
   FixedDiscretePolicy* vi_policy = value_iteration.getPolicy();
   PolicyEvaluation f_policy_evaluation(vi_policy, mdp, gamma);
   f_policy_evaluation.ComputeStateValuesFeatureExpectation(epsilon, -1);
   MDPWrapper wrapper(mdp);
   MDPWrapper* env = &wrapper;
   for (int s=0; s<mdp->getNStates(); s++) {
       int n_iter = 1000;
       int T = 1000;
       real ret = 0.0;
       real cumret = 0.0;
       for (int iter=0; iter<n_iter; ++iter) {
	 env->Reset(s);
	 real rt = env->getReward();
	 int st = env->getState();
	 vi_policy->Observe(rt, st);
	 real discount = 1.0;
	 for (int t=0; t<T; ++t) {
	   int at = vi_policy->SelectAction();
	   env->Act(at);
	   real rt = env->getReward();
	   int st = env->getState();
	   vi_policy->Observe(rt, st);
	   //printf("%f %d %d # rsa \n", rt, st, at);
	   ret += rt * discount;
	   cumret += rt;
	   discount *= gamma;
	 }
       }
       ret /= (real) n_iter;
       cumret /= (real) n_iter;
       printf ("%f %f %f ~ %f %f\n",
               value_iteration.getValue(s),
               policy_evaluation.getValue(s),
               f_policy_evaluation.getValue(s),
	       ret,
	       cumret);       
   }

   
   //mdp->ShowModel();
   if (1) {
       FILE* f = fopen("test.dot", "w");
       if (f) {
           fprintf (f, "digraph MDP {\n");
           fprintf (f, "ranksep=2; rankdir=LR; \n");
           mdp->dotModel(f);
           fprintf (f, "}\n");
           fclose(f);
       }
   }
#endif

#if 0
    // GRIDOWRLD CODE
    for (int iter=0; iter < 100; iter++) {
	printf ("ITER: %d\n", iter);
    for (uint y=0; y<height; y++) {
        for (uint x=0; x<width; x++) {
            int s= x + y*width;
            for (int a=0; a<n_actions; a++) {
                Q[a] = policy_evaluation.getValue(s,a);
            }
            //int a_max = ArgMax(Q);
            //real Q_max = Q[a_max];
            real V = policy_evaluation.getValue(s);
	    Gridworld::MapElement element = grid_world.whatIs(x, y);
	    if (x) printf ("&");
	    if (element == Gridworld::WALL) {
		printf ("\\bsqr ");
	    } else {
		printf ("%+.1f ", V);
	    }
        }
        printf ("\\\\\\hline\n");
    }


    for (uint y=0; y<height; y++) {
        for (uint x=0; x<width; x++) {
	    int s = grid_world.getState(x, y);
	    Gridworld::MapElement element = grid_world.whatIs(x, y);
		  if (x) printf ("&");

	    if (element == Gridworld::WALL) {
		//printf ("#");
		printf ("\\msqr ");
	    } else {
		for (int a=0; a<n_actions; a++) {
		    Q[a] = policy_evaluation.getValue(s,a);
		}
		int a_max = ArgMax(Q);
		switch (a_max) {
		case Gridworld::NORTH: printf ("$\\uparrow$ "); break;
		case Gridworld::SOUTH: printf ("$\\downarrow$ "); break;
		case Gridworld::EAST:  printf ("$\\rightarrow$ "); break;
		case Gridworld::WEST:  printf ("$\\leftarrow$ "); break;
		}
	    }
        }
        printf ("\\\\\\hline\n");
        //printf ("\n");
    }
    policy_evaluation.ComputeStateValues(0.001, 1);
    //policy_evaluation.ComputeStateActionValues(0.001, 1);
    }
    grid_world.Show();
#endif
}


#endif
