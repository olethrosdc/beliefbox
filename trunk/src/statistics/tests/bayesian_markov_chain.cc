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
#include "DiscreteHiddenMarkovModel.h"
#include "DiscreteHiddenMarkovModelEM.h"
#include "DiscreteHiddenMarkovModelPF.h"
#include "BayesianPredictiveStateRepresentation.h"
#include "BayesianMarkovChain.h"
#include "Random.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "EasyClock.h"
#include "Dirichlet.h"
#include "RandomDevice.h"
#include <ctime>
struct ErrorStatistics
{
    std::vector<real> loss;
    ErrorStatistics(int T) : loss(T)
    {}
};

void print_result(const char* fname, ErrorStatistics& error)
{
    FILE* f = fopen(fname, "w");
    if (!f) return;
    
    int T = error.loss.size();
    for (int t=0; t<T; ++t) {
        fprintf (f, " %f", error.loss[t]);
    }
    fprintf (f, "\n");
    fclose (f);
}

int main (int argc, char** argv)
{
    int n_observations = 4;
    int max_states = 4;
    int n_mc_states = 4;
    float prior = 1.0f;
    int T = 100;
    int n_iter = 100;
    //real stationarity = 0.99;
    
    setRandomSeed(time(NULL));

    if (argc > 1) {
        T = atoi(argv[1]);
    }

    if (argc > 2) {
        n_iter = atoi(argv[2]);
    }

    if (argc > 3) {
        n_observations = atoi(argv[3]);
    }

    if (argc > 4) {
        max_states = atoi(argv[4]);
    }

    if (argc > 5) {
        n_mc_states = atoi(argv[5]);
    }

    double oracle_time = 0;
    double bmc_time = 0;
    double bpsr_time = 0;
    double hmm_pf_time = 0;
    double hmm_is_pf_time = 0;
    double hmm_em_time = 0;
    double hmm_pf_mix_time = 0;

    double initial_time  = GetCPU();
    double elapsed_time = 0;

    ErrorStatistics oracle_error(T);
    ErrorStatistics bmc_error(T);
    ErrorStatistics bpsr_error(T);
    ErrorStatistics hmm_pf_error(T);
    ErrorStatistics hmm_is_pf_error(T);
    ErrorStatistics hmm_em_error(T);
    ErrorStatistics hmm_pf_mix_error(T);
    RandomDevice random_device(false);
    
    for (int iter=0; iter<n_iter; iter++) {
        //real stationarity = 0.9;
        real true_stationarity = 0.5 + 0.5*urandom();
        double remaining_time = (real) (n_iter - iter) * elapsed_time / (real) iter;
        printf ("# iter: %d, %.1f running, %.1f remaining\n", iter, elapsed_time, remaining_time);
        
        
        //logmsg ("Making Bayesian Markov chain\n");
        // our model for the chains
        bool dense = false;
        BayesianMarkovChain bmc(n_observations, 1+max_states, prior, dense);
        BayesianPredictiveStateRepresentation bpsr(n_observations, 1+max_states, prior, dense);

        //logmsg ("Making Markov chain\n");
        // the actual model that generates the data
        DiscreteHiddenMarkovModel* hmm = MakeRandomDiscreteHMM (n_mc_states,  n_observations,  true_stationarity, &random_device);
        DiscreteHiddenMarkovModelStateBelief oracle(hmm);

        real hmm_threshold = 0.5;
        real hmm_stationarity = 0.9;
        int hmm_particles = 128;
        DiscreteHiddenMarkovModelPF_ReplaceLowest hmm_pf(hmm_threshold, hmm_stationarity, n_mc_states, n_observations, hmm_particles);
        DiscreteHiddenMarkovModelPF_ISReplaceLowest hmm_is_pf(hmm_threshold, hmm_stationarity, n_mc_states, n_observations, hmm_particles);
        DiscreteHiddenMarkovModelEM hmm_em(n_mc_states, n_observations, hmm_stationarity, &random_device, 1);
        //DHMM_PF_Mixture<DiscreteHiddenMarkovModelPF> hmm_pf(hmm_threshold, hmm_stationarity, n_observations, hmm_particles, 2 * n_mc_states);
        DHMM_PF_Mixture<DiscreteHiddenMarkovModelPF_ReplaceLowest> hmm_pf_mix(hmm_threshold, hmm_stationarity, n_observations, hmm_particles, 2 * n_mc_states);

        //logmsg ("Observing chain outputs\n");
        oracle.Reset();
        bmc.Reset();
        bpsr.Reset();
        hmm->Reset();
        hmm_pf.Reset();
        hmm_is_pf.Reset();
        hmm_em.Reset();
        hmm_pf_mix.Reset();


        for (int t=0; t<T; ++t) {
            int observation = hmm->generate();

            int oracle_prediction = oracle.predict();
            if (oracle_prediction != observation) {
                oracle_error.loss[t] += 1.0;
            }

            int bmc_prediction = bmc.predict();
            if (bmc_prediction != observation) {
                bmc_error.loss[t] += 1.0;
            }

            int bpsr_prediction = bpsr.predict();
            if (bpsr_prediction != observation) {
                bpsr_error.loss[t] += 1.0;
            }

            int hmm_pf_prediction = hmm_pf.predict();
            if (hmm_pf_prediction != observation) {
                hmm_pf_error.loss[t] += 1.0;
            }

            int hmm_is_pf_prediction = hmm_is_pf.predict();
            if (hmm_is_pf_prediction != observation) {
                hmm_is_pf_error.loss[t] += 1.0;
            }

            int hmm_em_prediction = hmm_em.predict();
            if (hmm_em_prediction != observation) {
                hmm_em_error.loss[t] += 1.0;
            }

            int hmm_pf_mix_prediction = hmm_pf_mix.predict();
            if (hmm_pf_mix_prediction != observation) {
                hmm_pf_mix_error.loss[t] += 1.0;
            }
            
            double start_time, end_time;

            start_time = GetCPU();
            oracle.Observe(observation);
            end_time = GetCPU();
            oracle_time += end_time - start_time;

            start_time = end_time;
            bmc.ObserveNextState(observation);
            end_time = GetCPU();
            bmc_time += end_time - start_time;

            start_time = end_time;;
            bpsr.ObserveNextState(observation);                               
            end_time = GetCPU();
            bpsr_time += end_time - start_time;

            start_time = end_time;;
            hmm_pf.Observe(observation);                               
            end_time = GetCPU();
            hmm_pf_time += end_time - start_time;

            start_time = end_time;;
            hmm_is_pf.Observe(observation);                               
            end_time = GetCPU();
            hmm_is_pf_time += end_time - start_time;

            start_time = end_time;;
            hmm_em.Observe(observation);                               
            end_time = GetCPU();
            hmm_em_time += end_time - start_time;

            start_time = end_time;;
            hmm_pf_mix.Observe(observation);                               
            end_time = GetCPU();
            hmm_pf_mix_time += end_time - start_time;
        }

        double end_time = GetCPU();
        elapsed_time += end_time - initial_time;
        initial_time = end_time;

        delete hmm;
    }

    printf ("# Time -- Oracle: %f, HMM PF: %f, HMM IS PF: %f, HMM PF EX: %f, HMM PF MIX: %f, BHMC: %f, BPSR: %f\n", 
            oracle_time, hmm_pf_time, hmm_is_pf_time, hmm_em_time, hmm_pf_mix_time, bmc_time, bpsr_time);

    real inv_iter = 1.0 / (real) n_iter;
    for (int t=0; t<T; ++t) {
        hmm_pf_error.loss[t] *= inv_iter;
        hmm_is_pf_error.loss[t] *= inv_iter;
        hmm_em_error.loss[t] *= inv_iter;
        hmm_pf_mix_error.loss[t] *= inv_iter;
        oracle_error.loss[t] *= inv_iter;
        bmc_error.loss[t] *= inv_iter;
        bpsr_error.loss[t] *= inv_iter;
    }
    print_result("hmm_pf.error", hmm_pf_error);
    print_result("hmm_is_pf.error", hmm_is_pf_error);
    print_result("hmm_em.error", hmm_em_error);
    print_result("hmm_pf_mix.error", hmm_pf_mix_error);
    print_result("oracle.error", oracle_error);
    print_result("bmc.error", bmc_error);
    print_result("bpsr.error", bpsr_error);

}

#endif
