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
#include "DiscreteHiddenMarkovModelPF.h"
#include "BayesianPredictiveStateRepresentation.h"
#include "BayesianMarkovChain.h"
#include "Random.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "EasyClock.h"
#include "Dirichlet.h"
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

DiscreteHiddenMarkovModel* MakeRandomDiscreteHMM(int n_states, int n_observations, real stationarity)
{
    assert (n_states > 0);
    assert (n_observations > 0);
    assert (stationarity >= 0 && stationarity <= 1);

    Matrix Pr_S(n_states, n_states);
    Matrix Pr_X(n_states, n_observations);
    for (int i=0; i<n_states; ++i) {
        real sum = 0.0;
        for (int j=0; j<n_observations; ++j) {
            Pr_X(i,j) = 0.1*true_random(false);
            if (i==j) {
                Pr_X(i,j) += 1.0;
            }
            sum += Pr_X(i,j);
        }
        for (int j=0; j<n_observations; ++j) {
            Pr_X(i,j) /=  sum;
        }
        
    }
    Matrix S(n_states, n_states);
    for (int src=0; src<n_states; ++src) {
        Vector P(n_states);
        for (int i=0; i<n_states; ++i) {
            if (i<=src) {
                P[i] = exp(true_random(false));
            } else {
                P[i] = exp(10.0*true_random(false));
            }
        }
        P /= P.Sum();
        P *= (1 - stationarity);
        P[src] += stationarity;
        P /= P.Sum();
            //real sum = 0.0;
        for (int dst=0; dst<n_states; ++dst) {
            Pr_S(src,dst) = P[dst];
        }
    }

    return new DiscreteHiddenMarkovModel (Pr_S, Pr_X);
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
    double hmm_time = 0;
    double hmm_rep_time = 0;

    double initial_time  = GetCPU();
    double elapsed_time = 0;

    ErrorStatistics oracle_error(T);
    ErrorStatistics bmc_error(T);
    ErrorStatistics bpsr_error(T);
    ErrorStatistics hmm_error(T);
    ErrorStatistics hmm_rep_error(T);

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
        DiscreteHiddenMarkovModel* hmm = MakeRandomDiscreteHMM (n_mc_states,  n_observations,  true_stationarity);
        DiscreteHiddenMarkovModelStateBelief oracle(hmm);

        real hmm_threshold = 0.5;
        real hmm_stationarity = 0.9;
        int hmm_particles = 128;
        DiscreteHiddenMarkovModelPF_ReplaceLowestExact hmm_pf(hmm_threshold, hmm_stationarity, n_mc_states, n_observations, hmm_particles);
        //DHMM_PF_Mixture<DiscreteHiddenMarkovModelPF> hmm_pf(hmm_threshold, hmm_stationarity, n_observations, hmm_particles, 2 * n_mc_states);
        DHMM_PF_Mixture<DiscreteHiddenMarkovModelPF_ReplaceLowest> hmm_pf_rep(hmm_threshold, hmm_stationarity, n_observations, hmm_particles, 2 * n_mc_states);

        //logmsg ("Observing chain outputs\n");
        oracle.Reset();
        bmc.Reset();
        bpsr.Reset();
        hmm->Reset();
        

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

            int hmm_prediction = hmm_pf.predict();
            if (hmm_prediction != observation) {
                hmm_error.loss[t] += 1.0;
            }

            int hmm_rep_prediction = hmm_pf_rep.predict();
            if (hmm_rep_prediction != observation) {
                hmm_rep_error.loss[t] += 1.0;
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
            hmm_time += end_time - start_time;

            start_time = end_time;;
            hmm_pf_rep.Observe(observation);                               
            end_time = GetCPU();
            hmm_rep_time += end_time - start_time;
        }

        double end_time = GetCPU();
        elapsed_time += end_time - initial_time;
        initial_time = end_time;

        delete hmm;
    }

    printf ("# Time -- Oracle: %f, HMM: %f, HMM R: %f, BHMC: %f, BPSR: %f\n", 
            oracle_time, hmm_time, hmm_rep_time, bmc_time, bpsr_time);

    real inv_iter = 1.0 / (real) n_iter;
    for (int t=0; t<T; ++t) {
        hmm_error.loss[t] *= inv_iter;
        hmm_rep_error.loss[t] *= inv_iter;
        oracle_error.loss[t] *= inv_iter;
        bmc_error.loss[t] *= inv_iter;
        bpsr_error.loss[t] *= inv_iter;
    }
    print_result("hmm.error", hmm_error);
    print_result("hmm_rep.error", hmm_rep_error);
    print_result("oracle.error", oracle_error);
    print_result("bmc.error", bmc_error);
    print_result("bpsr.error", bpsr_error);

}

#endif
