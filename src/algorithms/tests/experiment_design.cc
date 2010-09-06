// -*- Mode: C++; -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/************************************************************************
 *                                                                      *
 * This program is free software; you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation; either version 2 of the License, or    *
 * (at your option) any later version.                                  *
 *                                                                      *
 ************************************************************************/

#ifdef MAKE_MAIN

#include "Random.h"
#include "RandomNumberFile.h"
#include "BetaDistribution.h"
#include <vector>
#include <cstdlib>

/* Consider the problem of online exploration for experimental design. 
   
   We are making a complete tree of beliefs.

*/

class Hypothesis;

class HypothesisEdge
{
public:
    Hypothesis* child;
    int action;
    real value;
    real probability;
};

/** Hypothesis.

 */
class Hypothesis
{
public:
    std::vector<HypothesisEdge*> edges;
    std::vector<BetaDistribution> belief;
    int n_actions;
    int t; ///< current tie
    int T; ///< horizon
    real Vmax; ///< value!
    std::vector<real> V; ///< Value of actions

    /// Create hypothesis form prior directly. Make children if necessary.
    Hypothesis(std::vector<BetaDistribution>& prior,
               int n_actions_,
               int t_,
               int T_,
               std::vector<real>& outcomes)

        : belief(prior), n_actions(n_actions_), t(t_), T(T_), Vmax(-INF), V(n_actions)
    {

        if (t < T) {
            MakeChildren(outcomes);
        } else {
            printf ("Ending at the prior! Whuh\n");
        }
    }
    
    /// Create hypothesis from prior and new observation. Make
    /// children if necessary.
    Hypothesis(std::vector<BetaDistribution>& prior,
               int n_actions_,
               int t_,
               int T_,
               std::vector<real>& outcomes,
               int action,
               real outcome)
        : belief(prior), n_actions(n_actions_), t(t_), T(T_), Vmax(-INF), V(n_actions)
    {
        prior[action].calculatePosterior(outcome);

        if (t < T) {
            MakeChildren(outcomes);
        } else {
            //printf ("Ending at action %d belief %f / %f\n", action, prior[action].alpha, prior[action].beta);
        }
    }

    /// How to make children!
    void MakeChildren(std::vector<real>& outcomes) {
        for (int a=0; a<n_actions; ++a) {
            for (uint i=0; i!=outcomes.size(); ++i) {
                real u = belief[a].getMean();
                real x = outcomes[i];
                real p = u * x + (1 - u) * x;
                HypothesisEdge* edge = new HypothesisEdge;
                edge->child = new Hypothesis(belief,
                                             n_actions,
                                             t + 1,
                                             T,
                                             outcomes,
                                             a,
                                             i);
                edge->action = a;
                edge->probability = p;
                edge->value = 0;
                //printf("New edge: a=%d, x=%f, p=%f\n", a, outcomes[i], p);
                edges.push_back(edge);
            }
        }
    }

    /// Get the value, looking at children if necessary.
    real Value()
    {
        if (t == T) {
            return ExpectedValue();
        } 


        for (int a=0; a<n_actions; ++a) {
            V[a] = 0;
        }

        for (uint i=0; i<edges.size(); ++i) {
            edges[i]->value = edges[i]->child->Value();
            V[edges[i]->action] += edges[i]->value * edges[i]->probability;
#if 0
            printf(" - edge %d = %f, wp %f\n",
                   edges[i]->action,
                   edges[i]->value,
                   edges[i]->probability);
#endif

        }
        Vmax = Max(V);
        //printf ("=> %f\n", Vmax);
        return Vmax;
    }

    /// Get the value of making a decision right now!
    real ExpectedValue()
    {
        Vmax = -INF;
        for (int a=0; a<n_actions; ++a) {
            real p = belief[a].alpha / (belief[a].alpha + belief[a].beta);
            V[a] = p * 1.0 + (-10) * p;
            if (V[a] > Vmax) {
                Vmax = V[a];
            }
        }
        //printf ("EV: %f\n", Vmax);
        return Vmax;
    }

    void ShowValues()
    {
        for (int a=0; a<n_actions; ++a) {
            V[a] = 0;
        }

        for (uint i=0; i<edges.size(); ++i) {
            edges[i]->value = edges[i]->child->Value();
            V[edges[i]->action] += edges[i]->value * edges[i]->probability;
        }

        for (int a=0; a<n_actions; ++a) {
            printf ("V(%d) = %f\n", a, V[a]);
        }

    }
};


class ExperimentDesign
{
public:
    Hypothesis* root;
    std::vector<real>& outcomes;
    int n_actions;
    std::vector<BetaDistribution> prior;
    int T;
    ExperimentDesign(std::vector<real>& outcomes_,
                     int n_actions_,
                     int T_,
                     real alpha,
                     real beta)
        : outcomes(outcomes_),
          n_actions(n_actions_),
          prior(n_actions),
          T(T_)
    {
        for (int a=0; a<n_actions; ++a)  {
            prior[a].alpha = 0.5 + a;
            prior[a].beta = 0.5;
        }
        root = new Hypothesis(prior,
                              n_actions,
                              1,
                              T,
                              outcomes);
        
    }

    real Value()
    {
        return root->Value();
    }

    void ShowValues()
    {
        root->ShowValues();
    }
    
};

int main(int argc, char** argv)
{
    int n_actions = atoi(argv[1]);
    int T = atoi(argv[2]) ;

    if (argc != 3) {
        fprintf (stderr, "Usage: experiment_design n_actions T\n");
        return -1;
    }

    if (n_actions < 1) {
        fprintf (stderr, "n_actions > 0\n");
        return -1;
    }

    if (T < 1) {
        fprintf (stderr, "T > 0\n");
        return -1;
    }

    std::vector<real> outcomes(2);
    outcomes[0] = 0;
    outcomes[1] = 1;

    ExperimentDesign experiment_design(outcomes, n_actions, T, 1.0, 1.0);
    printf ("Value : %f\n", experiment_design.Value());
    experiment_design.ShowValues();
    return 0;
}

#endif
