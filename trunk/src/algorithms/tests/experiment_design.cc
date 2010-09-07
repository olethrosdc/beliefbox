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

#undef DEBUG_EXPERIMENT

#ifdef MAKE_MAIN

#include "Random.h"
#include "RandomNumberFile.h"
#include "BetaDistribution.h"
#include "RandomDevice.h"
#include "RandomNumberGenerator.h"
#include <vector>
#include <cstdlib>

#define CORRECT_VALUE 1
#define INCORRECT_VALUE 0

/* Consider the problem of online exploration for experimental design. 
   
   We are making a complete tree of beliefs.

*/

class Hypothesis;
class HypothesisEdge {
public:
    Hypothesis* child;
    int action;
    real value;
    real probability;
    ~HypothesisEdge();
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
#ifdef DEBUG_EXPERIMENT
        for (uint i=0; i<belief.size(); ++i) {
            printf ("B %d %f\n", i, belief[i].getMean());
        }
#endif

        if (t < T) {
            MakeChildren(outcomes);
        } else {
            //printf ("Ending at the prior! Whuh\n");
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

    ~Hypothesis()
    {
        for (uint i=0; i<edges.size(); ++i) {
            delete edges[i];
        }
    }
    /// How to make children!
    void MakeChildren(std::vector<real>& outcomes) {
        for (int a=0; a<n_actions; ++a) {
            real u = belief[a].getMean();
            //printf ("u[%d] = %f\n", a, u);
            for (uint i=0; i!=outcomes.size(); ++i) {
                real x = outcomes[i];
                real p = u * x + (1 - u) * (1 - x);
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
                //printf("New edge: t = %d, a=%d, x=%f, p=%f\n", t, a, outcomes[i], p);
                edges.push_back(edge);
            }
        }
    }

    /// Get the value, looking at children if necessary.
    real Value()
    {
        if (t == T) {
            Vmax = ExpectedValue();
            return Vmax;
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
            V[a] = p * CORRECT_VALUE + INCORRECT_VALUE * (1 - p);
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
            for (int k=0; k<t; ++k) {
                printf ("-");
            }
            printf ("[%d (%d): %f w.p. %f]\n", 
                    i,
                    edges[i]->action,
                    edges[i]->value,
                    edges[i]->probability);
            edges[i]->child->ShowValues();
        }
        for (int a=0; a<n_actions; ++a) {
            for (int k=0; k<t; ++k) {
                printf ("=");
            }
            printf ("V(%d) = %f\n", a, V[a]);
        }

    }
    int SelectAction()
    {
        return ArgMax(V);
    }
};

HypothesisEdge::~HypothesisEdge()
{
    delete child;
}

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
                     std::vector<BetaDistribution>& prior_,
                     int T_)
        : outcomes(outcomes_),
          n_actions(n_actions_),
          prior(prior_),
          T(T_)
    {
        root = new Hypothesis(prior,
                              n_actions,
                              0,
                              T,
                              outcomes);
        
    }
    ~ExperimentDesign()
    {
        delete root;
    }
    real Value()
    {
        return root->Value();
    }

    void ShowValues()
    {
        root->ShowValues();
    }
    int SelectAction()
    {
        return root->SelectAction();
    }
};
void Experiment(int n_actions, int T);
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
    Experiment (n_actions, T);
    return 0;
}

void Experiment(int n_actions, int T)
{
    std::vector<real> outcomes(2);
    outcomes[0] = 0;
    outcomes[1] = 1;
    
    int n_iter = 1000;
    real randomised_results = 0;
    real lookahead_results = 0;
    for (int k=0; k<n_iter; ++k) {
        std::vector<real> success_probabilities(2);
        for (int i=0; i<n_actions; ++i) {
            success_probabilities[i] = urandom();//0.25 + 0.75 * ((real) i) / ((real) n_actions);
#ifdef DEBUG_EXPERIMENT
            printf ("P[%d] = %f\n", i, success_probabilities[i]);
#endif
        }

        
        // Clever trials
        {
            std::vector<BetaDistribution> prior(n_actions);
            for (int a=0; a<n_actions; ++a)  {
                prior[a].alpha = 1.0;
                prior[a].beta = 1.0;
            }
            for (int t=0; t < T; ++t) {
                ExperimentDesign experiment_design(outcomes, n_actions, prior, T - t);                
                int action = experiment_design.SelectAction();
                real outcome;
                if (urandom() < success_probabilities[action]) {
                    outcome = outcomes[1];
                } else {
                    outcome = outcomes[0];
                }
#ifdef DEBUG_EXPERIMENT
                printf ("Value : %f\n", experiment_design.Value());
                experiment_design.ShowValues();
                printf("Action : %d\n", action);
                printf ("Outcome: %f\n", outcome);
#endif
                prior[action].alpha += outcome;
                prior[action].beta += (1 - outcome);
            }
            std::vector<real> V(n_actions);
            for (int a=0; a<n_actions; ++a) {
                real p = prior[a].alpha / (prior[a].alpha + prior[a].beta);
                V[a] = p * CORRECT_VALUE + (1 - p) * INCORRECT_VALUE;
            }
            int a_max = ArgMax(V);
#ifdef DEBUG_EXPERIMENT
            printf ("Selected Treatment : %d\n", a_max);
#endif
            real true_value = success_probabilities[a_max] * CORRECT_VALUE
                + (1 - success_probabilities[a_max]) * INCORRECT_VALUE;
            lookahead_results += true_value;
        }

        // Randomised trials
        {
            std::vector<BetaDistribution> prior(n_actions);
            for (int a=0; a<n_actions; ++a)  {
                prior[a].alpha = 1.0;
                prior[a].beta = 1.0;
            }
            RandomDevice rng(false);
            for (int t=0; t < T; ++t) {
                int action = rng.discrete_uniform(n_actions);
                //int action = t % n_actions;
                real outcome;
                if (urandom() < success_probabilities[action]) {
                    outcome = outcomes[1];
                } else {
                    outcome = outcomes[0];
                }
                //printf ("Outcome: %f\n", outcome);
                prior[action].alpha += outcome;
                prior[action].beta += (1 - outcome);
            }
            std::vector<real> V(n_actions);
            for (int a=0; a<n_actions; ++a) {
                real p = prior[a].alpha / (prior[a].alpha + prior[a].beta);
                V[a] = p * CORRECT_VALUE + (1 - p) * INCORRECT_VALUE;
            }
            int a_max = ArgMax(V);
            //printf ("Selected Treatment : %d\n", a_max);
            real true_value = success_probabilities[a_max] * CORRECT_VALUE
                + (1 - success_probabilities[a_max]) * INCORRECT_VALUE;
            randomised_results += true_value;
        }
    }
    printf ("%f %f #VAL\n",
            lookahead_results / (real) n_iter,
            randomised_results / (real) n_iter);
        
}

#endif
