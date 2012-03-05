/* -*- Mode: C++; -*- */
/* VER: $Id: test_vpi_evaluation.c,v 1.12 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $*/
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
#include "UCB.h"

#include "ActionValueEstimate.h"
#include "PFActionValueEstimate.h"
#include "EasyClock.h"

void Evaluate(char* fname, DiscreteBanditPolicy& policy, int n_actions);

// if the result is worse but it seems that it could be better then we have to
// calculate that.
int main (int argc, char** argv)
{
    fprintf (stderr, "test_pf_evaluation $Id: test_vpi_evaluation.c,v 1.12 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $\n");
    if (argc!=4) {
        fprintf (stderr, "Args: n-actions n_samples n_samples2\n");
        exit(-1);
    }
        
    int n_actions = atoi(argv[1]);
    int n_samples = atoi(argv[2]);
    int n_samples2 = atoi(argv[3]);
    fprintf (stderr, "n_actions: %d n_samples: %d n_samples2: %d\n", n_actions, n_samples, n_samples2);
    //real alpha = 0.1;
    real alpha = 0.1;
    real epsilon  = 0.01;
    real gamma = 0.999;
    PopulationEstimate population_estimate(n_actions, n_samples, alpha);
    PFActionValueEstimate pf_estimate(n_actions, n_samples);
    BernoulliEstimate bernoulli_estimate (n_actions, n_samples, 1.0f, 1.0f);

    srand(time(NULL));
    srand48(time(NULL));

    PopulationVPIPolicy vpi_policy0(n_actions, &population_estimate, 0.0);
    PopulationVPIPolicy vpi_policy5(n_actions, &population_estimate, 0.5);
    PopulationVPIPolicy vpi_policy9(n_actions, &population_estimate, 0.9);
    PopulationVPIPolicy vpi_policy99(n_actions, &population_estimate, 0.99);
    PopulationVPIPolicy vpi_policy999(n_actions, &population_estimate, 0.999);
    VPIPolicy bayes_vpi_policy0(n_actions, &bernoulli_estimate,0.0,n_samples2);
    VPIPolicy bayes_vpi_policy5(n_actions, &bernoulli_estimate,0.5,n_samples2);
    VPIPolicy bayes_vpi_policy9(n_actions, &bernoulli_estimate,0.9,n_samples2);
    VPIPolicy bayes_vpi_policy99(n_actions, &bernoulli_estimate,0.99,n_samples2);
    VPIPolicy bayes_vpi_policy999(n_actions, &bernoulli_estimate,0.999,n_samples2);

    PFVPIPolicy pf_vpi_policy(n_actions, &pf_estimate, gamma);

    PopulationSamplePolicy bayes_sample(n_actions, &bernoulli_estimate, 1.0);
    PopulationOptimalPolicy bayes_policy5(n_actions, &bernoulli_estimate, 0.5, n_samples2);
    PopulationOptimalPolicy bayes_policy9(n_actions, &bernoulli_estimate, 0.9, n_samples2);
    PopulationOptimalPolicy bayes_policy99(n_actions, &bernoulli_estimate, 0.99, n_samples2);
    PopulationOptimalPolicy bayes_policy999(n_actions, &bernoulli_estimate, 0.999, n_samples2);
    PopulationOptimalPolicy bayes_policy9999(n_actions, &bernoulli_estimate, 0.9999, n_samples2);
    OptimalInfinitePolicy optimal_infinite(n_actions);
    UCB1Policy ucb1(n_actions);
    UCB2Policy ucb2a1(0.1, n_actions);
    UCB2Policy ucb2a25(0.25, n_actions);
    UCB2Policy ucb2a5(0.5, n_actions);
    UCB2Policy ucb2a75(0.75, n_actions);
    UCB2Policy ucb2a9(0.9, n_actions);
    UCB2Policy ucb2a99(0.99, n_actions);

    //Evaluate ("out_opt_inf", optimal_infinite, n_actions);
    Evaluate ("out_bayes_sample", bayes_sample, n_actions);
    Evaluate ("out_ucb2a99", ucb2a99, n_actions);
    Evaluate ("out_ucb2a9", ucb2a9, n_actions);
    Evaluate ("out_ucb2a75", ucb2a75, n_actions);
    Evaluate ("out_ucb1", ucb1, n_actions);
    Evaluate ("out_ucb2a1", ucb2a1, n_actions);
    Evaluate ("out_ucb2a25", ucb2a25, n_actions);
    Evaluate ("out_ucb2a5", ucb2a5, n_actions);

    Evaluate ("out_bay999_cor", bayes_policy999, n_actions);
    Evaluate ("out_bay5_cor", bayes_policy5, n_actions);
    Evaluate ("out_bay9_cor", bayes_policy9, n_actions);
    Evaluate ("out_bay99_cor", bayes_policy99, n_actions);

    Evaluate ("out_bay9999_cor", bayes_policy9999, n_actions);
    Evaluate ("out_bay_vpi999", bayes_vpi_policy999, n_actions);
    Evaluate ("out_pf_vpi999", pf_vpi_policy, n_actions);

    return 0;
}


void Evaluate(char* fname, DiscreteBanditPolicy& policy, int n_actions)
{
    int episode_length = 100000;
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


    double elapsed_time=0.0;
    double start_time = GetCPU();
    for (int episode = 0; episode < n_episodes; episode++) {
        policy.Reset();
        std::vector<real> Er(n_actions);
        for (int i=0; i<n_actions; i++) {
            Er[i] = true_random(false);
        }
        for (int t=0; t<episode_length; t++) {
            int a = policy.SelectAction();
            real r = 0.0;
            if (urandom() < Er[a]) {
                r = 1.0;
            }
            stats[t] += r;
            policy.Observe(a,r);
        }


        elapsed_time = GetCPU() - start_time;
        double time_per_episode = (double) elapsed_time / (double) (episode+1);
        double total_time = time_per_episode * (double) n_episodes;
        double remaining_time = total_time - (double) elapsed_time;
        double percent_complete = 100.0f * elapsed_time/total_time;
        if (std::isnan(percent_complete) || percent_complete<0) {
            percent_complete = 100.0f * (double) (episode+1) / (double) n_episodes;
        }
        printf ("Completion %d%%. Time: %d elapsed, %d total, %d remaining.\n",
                (int) percent_complete,
                (int) elapsed_time,
                (int) total_time,
                (int) remaining_time);
        
    }
    
    printf ("%f seconds passed\n", elapsed_time);
    fprintf (f, "# %f seconds passed\n",  elapsed_time);
    for (int t=0; t<episode_length; t++) {
        stats[t] /= (double) n_episodes;
        fprintf (f, "%f\n", stats[t]);
    }
    fclose (f);
}

#endif
