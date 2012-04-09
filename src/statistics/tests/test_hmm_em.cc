/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "DiscreteHiddenMarkovModel.h"
#include "DiscreteHiddenMarkovModelEM.h"
#include "ExpectationMaximisation.h"
#include "Matrix.h"
#include "Dirichlet.h"
#include "Random.h"
#include "RandomDevice.h"
#include "CumulativeStats.h"

template <typename T>
void PredictAndAdapt(T& model, int x, int t, CumulativeStats& stats)
{
    int prediction = ArgMax(model.getPrediction());
    real loss = 1;
    if (prediction == x) {
        loss = 0;
    }
    stats.SetValue(t, loss);
    model.Observe(x);
}

void TestBelief (DiscreteHiddenMarkovModel* hmm,
                 int T,
                 real threshold,
                 real stationarity,
                 int n_em_iter,
                 CumulativeStats& oracle_state_stats,
                 CumulativeStats& em_state_stats)
{
    RandomDevice rng(false);
    DiscreteHiddenMarkovModelStateBelief hmm_belief_state(hmm);
    DiscreteHiddenMarkovModel* estimated_hmm_ptr = MakeRandomDiscreteHMM(hmm->getNStates(), hmm->getNObservations(), stationarity, &rng);
    DiscreteHiddenMarkovModel& estimated_hmm = *estimated_hmm_ptr;

    ExpectationMaximisation<DiscreteHiddenMarkovModel, int> EM_algo(estimated_hmm);
    //DiscreteHiddenMarkovModelEM EM_hmm(hmm->getNStates(), hmm->getNObservations(), stationarity, &rng, n_em_iter);
    std::vector<int> s(T);
    for (int t=0; t<T; ++t) {
        // generate next observation 
        int x = hmm->generate();
        s[t] = hmm->getCurrentState();

        // oracle
        int oracle_x = hmm_belief_state.predict();
        real oracle_loss = 0;
        if (oracle_x != x) {
            oracle_loss = 1;
        }
        real oracle_accuracy = hmm_belief_state.getPrediction()(x);

        // EM HMM
        int em_hmm_x = EM_hmm.predict();
        real em_hmm_loss = 0;
        if (em_hmm_x != x) {
            em_hmm_loss = 1;
        }
        real em_hmm_accuracy = EM_hmm.getPrediction()(x);

        // save stats
        oracle_state_stats.SetValue(t, oracle_accuracy);
        em_state_stats.SetValue(t, em_hmm_accuracy);
        
        //EM_algo.Observe(x);
        //if (t > 0) {
        //    EM_algo.Iterate(1);
        //}
        hmm_belief_state.Observe(x);
        EM_hmm.Observe(x);
        //oracle_state_stats.SetValue(t, hmm_belief_state.getBelief()[s[t]]);
    }

        //    EM_algo.Iterate(n_em_iter);
    //Matrix& belief = estimated_hmm.getBelief();
    //for (int t=0; t<T; ++t) {
    //        em_state_stats.SetValue(t, belief(t, s[t]));
        //printf ("%d %f\n", t, belief(t, s[t]));
    //}
    delete estimated_hmm_ptr;
}


int main(int argc, char** argv)
{
    if (argc != 7) {
        fprintf(stderr, "Usage: test_hmm_em n_states n_observations stationarity  EM_iter T n_iter\n");
        return -1;
    }
    int n_states = atoi(argv[1]);
    if (n_states <= 0) {
        fprintf (stderr, "Invalid number of states %d\n", n_states);
    }

    int n_observations = atoi(argv[2]);
    if (n_observations <= 0) {
        fprintf (stderr, "Invalid number of states %d\n", n_observations);
    }

    real stationarity = atof(argv[3]);
    if (stationarity < 0 || stationarity > 1) {
        fprintf (stderr, "Invalid stationarity %f\n", stationarity);
    }
    int n_em_iter = atoi(argv[4]);
    if (n_em_iter <= 0) {
        fprintf (stderr, "Invalid n_em_iter %d\n", n_em_iter);
    }

    int T = atoi(argv[5]);
    if (T <= 0) {
        fprintf (stderr, "Invalid T %d\n", T);
    }

    int n_iter = atoi(argv[6]);
    if (n_iter <= 0) {
        fprintf (stderr, "Invalid n_iter %d\n", n_iter);
    }
    

    real threshold = 0.5; // threshold for the prior in the estimated HMMs.

    Vector x(10);
    Vector y = exp(x);


    CumulativeStats oracle_state_stats(T, n_iter);
    CumulativeStats em_state_stats(T, n_iter);

    RandomDevice random_device(false);

    for (int i=0; i<n_iter; ++i) {
        real true_stationarity = 0.5 + 0.5 * urandom();
        oracle_state_stats.SetSequence(i);
        em_state_stats.SetSequence(i);
        //fprintf (stderr, "Iter: %d / %d\n", i + 1, n_iter);
        DiscreteHiddenMarkovModel* hmm = MakeRandomDiscreteHMM(n_states,  n_observations, true_stationarity, &random_device);
        TestBelief(hmm, T, threshold, stationarity, n_em_iter, oracle_state_stats, em_state_stats);
        delete hmm;
    }
    
    real percentile = 0.1;
    oracle_state_stats.Sort();
    Vector oracle_state_mean = oracle_state_stats.Mean();
    Vector oracle_state_top = oracle_state_stats.TopPercentile(percentile);
    Vector oracle_state_bottom = oracle_state_stats.BottomPercentile(percentile);

    em_state_stats.Sort();
    Vector em_state_mean = em_state_stats.Mean();
    Vector em_state_top = em_state_stats.TopPercentile(percentile);
    Vector em_state_bottom = em_state_stats.BottomPercentile(percentile);

    
    for (int t=0; t<T; ++t) {
        printf ("%f %f %f %f %f %f\n",
                em_state_bottom[t], em_state_mean[t], em_state_top[t],
                oracle_state_bottom[t], oracle_state_mean[t], oracle_state_top[t]);
    }
    return 0;
}

#endif
