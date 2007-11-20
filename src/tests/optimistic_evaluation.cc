/* -*- Mode: C++; -*- */
/* VER: $Id: test_optimistic_evaluation.c,v 1.3 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "DiscreteBanditPolicy.h"
#include "ActionValueEstimate.h"


void Evaluate(char* fname, Policy& policy, int n_actions);

// if the result is worse but it seems that it could be better then we have to
// calculate that.
int main (void)
{
 
    int n_actions = 16;
    int n_samples = 16;
    //real alpha = 0.1;
    real alpha = 0.1;
    real epsilon  = 0.01;
    real gamma = 0.9;
    PopulationEstimate population_estimate(n_actions, n_samples, alpha);
    PopulationEstimate population_estimate2(n_actions, n_samples, alpha);
    PointEstimate point_estimate (n_actions, alpha);
    PointEstimate optimistic_point_estimate_a001 (n_actions, 0.001f, 1.0f, 0.0f);
    PointEstimate optimistic_point_estimate_a01 (n_actions, 0.01f, 1.0f, 0.0f);
    PointEstimate optimistic_point_estimate_a1 (n_actions, 0.1f, 1.0f, 0.0f);

    srand(time(NULL));
    srand48(time(NULL));


    PopulationOptimalPolicy optimal_policy1_1(n_actions, &population_estimate, 0.1, 1);
    PopulationOptimalPolicy optimal_policy1_2(n_actions, &population_estimate, 0.1, 2);
    PopulationOptimalPolicy optimal_policy1_4(n_actions, &population_estimate, 0.1, 4);
    PopulationOptimalPolicy optimal_policy1_8(n_actions, &population_estimate, 0.1, 8);
    PopulationOptimalPolicy optimal_policy1_16(n_actions, &population_estimate, 0.1, 16);
    PopulationOptimalPolicy optimal_policy1_32(n_actions, &population_estimate, 0.1, 32);

    PopulationOptimalPolicy optimal_policy5_1(n_actions, &population_estimate, 0.5, 1);
    PopulationOptimalPolicy optimal_policy5_2(n_actions, &population_estimate, 0.5, 2);
    PopulationOptimalPolicy optimal_policy5_4(n_actions, &population_estimate, 0.5, 4);
    PopulationOptimalPolicy optimal_policy5_8(n_actions, &population_estimate, 0.5, 8);
    PopulationOptimalPolicy optimal_policy5_16(n_actions, &population_estimate, 0.5, 16);
    PopulationOptimalPolicy optimal_policy5_32(n_actions, &population_estimate, 0.5, 32);
    PopulationOptimalPolicy optimal_policy9_1(n_actions, &population_estimate, 0.9, 1);
    PopulationOptimalPolicy optimal_policy9_2(n_actions, &population_estimate, 0.9, 2);
    PopulationOptimalPolicy optimal_policy9_4(n_actions, &population_estimate, 0.9, 4);
    PopulationOptimalPolicy optimal_policy9_8(n_actions, &population_estimate, 0.9, 8);
    PopulationOptimalPolicy optimal_policy9_16(n_actions, &population_estimate, 0.9, 16);
    PopulationOptimalPolicy optimal_policy9_32(n_actions, &population_estimate, 0.9, 32);
    PopulationOptimalPolicy optimal_policy99_1(n_actions, &population_estimate, 0.99, 1);
    PopulationOptimalPolicy optimal_policy99_2(n_actions, &population_estimate, 0.99, 2);
    PopulationOptimalPolicy optimal_policy99_4(n_actions, &population_estimate, 0.99, 4);
    PopulationOptimalPolicy optimal_policy99_8(n_actions, &population_estimate, 0.99, 8);
    PopulationOptimalPolicy optimal_policy99_16(n_actions, &population_estimate, 0.99, 16);
    PopulationOptimalPolicy optimal_policy99_32(n_actions, &population_estimate, 0.99, 32);
    PopulationOptimalPolicy optimal_policy999_1(n_actions, &population_estimate, 0.999, 1);
    PopulationOptimalPolicy optimal_policy999_2(n_actions, &population_estimate, 0.999, 2);
    PopulationOptimalPolicy optimal_policy999_4(n_actions, &population_estimate, 0.999, 4);
    PopulationOptimalPolicy optimal_policy999_8(n_actions, &population_estimate, 0.999, 8);
    PopulationOptimalPolicy optimal_policy999_16(n_actions, &population_estimate, 0.999, 16);
    PopulationOptimalPolicy optimal_policy999_32(n_actions, &population_estimate, 0.999, 32);

    PopulationSamplePolicy sample_policy(n_actions, &population_estimate2, gamma);
    PopulationVPIPolicy vpi_policy(n_actions, &population_estimate2, gamma);

    EpsilonGreedyPolicy greedy_policy(n_actions, 0.0, &point_estimate);
    EpsilonGreedyPolicy egreedy_policy001(n_actions, 0.001, &point_estimate);
    EpsilonGreedyPolicy egreedy_policy01(n_actions, 0.01, &point_estimate);
    EpsilonGreedyPolicy egreedy_policy1(n_actions, 0.1, &point_estimate);

    EpsilonGreedyPolicy a001greedy_policy(n_actions, 0.0, &optimistic_point_estimate_a001);
    EpsilonGreedyPolicy a001egreedy_policy001(n_actions, 0.001, &optimistic_point_estimate_a001);
    EpsilonGreedyPolicy a001egreedy_policy01(n_actions, 0.01, &optimistic_point_estimate_a001);
    EpsilonGreedyPolicy a001egreedy_policy1(n_actions, 0.1, &optimistic_point_estimate_a001);

    EpsilonGreedyPolicy a01greedy_policy(n_actions, 0.0, &optimistic_point_estimate_a01);
    EpsilonGreedyPolicy a01egreedy_policy001(n_actions, 0.001, &optimistic_point_estimate_a01);
    EpsilonGreedyPolicy a01egreedy_policy01(n_actions, 0.01, &optimistic_point_estimate_a01);
    EpsilonGreedyPolicy a01egreedy_policy1(n_actions, 0.1, &optimistic_point_estimate_a01);

    EpsilonGreedyPolicy a1greedy_policy(n_actions, 0.0, &optimistic_point_estimate_a1);
    EpsilonGreedyPolicy a1egreedy_policy001(n_actions, 0.001, &optimistic_point_estimate_a1);
    EpsilonGreedyPolicy a1egreedy_policy01(n_actions, 0.01, &optimistic_point_estimate_a1);
    EpsilonGreedyPolicy a1egreedy_policy1(n_actions, 0.1, &optimistic_point_estimate_a1);




    NaiveE3Policy naive_e3_policy01(n_actions, epsilon, 0.01);
    NaiveE3Policy naive_e3_policy1(n_actions, epsilon, 0.1);
    NaiveE3Policy naive_e3_policy25(n_actions, epsilon, 0.25);
    NaiveE3Policy naive_e3_policy5(n_actions, epsilon, 0.5);
    NaiveE3Policy naive_e3_policy75(n_actions, epsilon, 0.75);
    NaiveE3Policy naive_e3_policy9(n_actions, epsilon, 0.9);
    NaiveE3Policy naive_e3_policy99(n_actions, epsilon, 0.99);
        

    Evaluate ("a1o_grd", a1greedy_policy, n_actions);
    Evaluate ("a1o_egr001", a1egreedy_policy001, n_actions);
    Evaluate ("a1o_egr01", a1egreedy_policy01, n_actions);
    Evaluate ("a1o_egr1", a1egreedy_policy1, n_actions);

    Evaluate ("a01o_grd", a01greedy_policy, n_actions);
    Evaluate ("a01o_egr001", a01egreedy_policy001, n_actions);
    Evaluate ("a01o_egr01", a01egreedy_policy01, n_actions);
    Evaluate ("a01o_egr1", a01egreedy_policy1, n_actions);

    Evaluate ("a001o_grd", a001greedy_policy, n_actions);
    Evaluate ("a001o_egr001", a001egreedy_policy001, n_actions);
    Evaluate ("a001o_egr01", a001egreedy_policy01, n_actions);
    Evaluate ("a001o_egr1", a001egreedy_policy1, n_actions);

#if 0
    Evaluate ("a1_grd", greedy_policy, n_actions);
    Evaluate ("a1_egr001", egreedy_policy001, n_actions);
    Evaluate ("a1_egr01", egreedy_policy01, n_actions);
    Evaluate ("a1_egr1", egreedy_policy1, n_actions);


    Evaluate ("out_sam", sample_policy, n_actions);
    Evaluate ("out_vpi", vpi_policy, n_actions);


    Evaluate ("out_opt1_1", optimal_policy1_1, n_actions);
    Evaluate ("out_opt1_2", optimal_policy1_2, n_actions);
    Evaluate ("out_opt1_4", optimal_policy1_4, n_actions);
    //Evaluate ("out_opt1_8", optimal_policy1_8, n_actions);
    //Evaluate ("out_opt1_16", optimal_policy1_16, n_actions);
    //Evaluate ("out_opt1_32", optimal_policy1_32, n_actions);

    Evaluate ("out_opt5_1", optimal_policy5_1, n_actions);
    Evaluate ("out_opt5_2", optimal_policy5_2, n_actions);
    Evaluate ("out_opt5_4", optimal_policy5_4, n_actions);
    Evaluate ("out_opt5_8", optimal_policy5_8, n_actions);
    //Evaluate ("out_opt5_16", optimal_policy5_16, n_actions);
    //Evaluate ("out_opt5_32", optimal_policy5_32, n_actions);

    Evaluate ("out_opt9_1", optimal_policy9_1, n_actions);
    Evaluate ("out_opt9_2", optimal_policy9_2, n_actions);
    Evaluate ("out_opt9_4", optimal_policy9_4, n_actions);
    //Evaluate ("out_opt9_8", optimal_policy9_8, n_actions);
    //Evaluate ("out_opt9_16", optimal_policy9_16, n_actions);
    //Evaluate ("out_opt9_32", optimal_policy9_32, n_actions);

    Evaluate ("out_opt99_1", optimal_policy99_1, n_actions);
    Evaluate ("out_opt99_2", optimal_policy99_2, n_actions);
    Evaluate ("out_opt99_4", optimal_policy99_4, n_actions);
    //Evaluate ("out_opt99_8", optimal_policy99_8, n_actions);
    //Evaluate ("out_opt99_16", optimal_policy99_16, n_actions);
    //Evaluate ("out_opt99_32", optimal_policy99_32, n_actions);

    Evaluate ("out_e3n01", naive_e3_policy01, n_actions);
    Evaluate ("out_e3n1", naive_e3_policy1, n_actions);
    Evaluate ("out_e3n25", naive_e3_policy25, n_actions);
    Evaluate ("out_e3n5", naive_e3_policy5, n_actions);
    Evaluate ("out_e3n75", naive_e3_policy75, n_actions);
    Evaluate ("out_e3n9", naive_e3_policy9, n_actions);
    Evaluate ("out_e3n99", naive_e3_policy99, n_actions);
#endif
    return 0;
}

void Evaluate(char* fname, Policy& policy, int n_actions)
{
    int episode_length = 10000;
    int n_episodes = 1000;

    printf ("Writing to %s\n", fname);
    fflush(stdout);

    FILE* f = fopen (fname, "w");
    if (!f) {
        fprintf (stderr, "Could not open %s for writing\n", fname);
        return;
    }


    std::vector<real> stats(episode_length);
    for (int t=0; t<episode_length; t++) {
        stats[t] = 0.0f;
    }

    time_t start_time = time(NULL);

    for (int episode = 0; episode < n_episodes; episode++) {
        policy.Reset();
        std::vector<real> Er(n_actions);
        for (int i=0; i<n_actions; i++) {
            Er[i] = drand48();
        }
        for (int t=0; t<episode_length; t++) {
            int a = policy.SelectAction();
            real r = 0.0;
            if (drand48() < Er[a]) {
                r = 1.0;
            }
            stats[t] += r;
            policy.Observe(a,r);
        }
        
    }
    
    time_t elapsed_time = time(NULL) - start_time;
    printf ("%d seconds passed\n", (int) elapsed_time);
    fprintf (f, "# %d seconds passed\n", (int) elapsed_time);
    for (int t=0; t<episode_length; t++) {
        stats[t] /= (real) n_episodes;
        fprintf (f, "%f\n", stats[t]);
    }
    fclose (f);
}

#endif
