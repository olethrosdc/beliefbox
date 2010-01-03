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
#include "BVMM.h"
#include "BayesianHierarchicalMarkovChain.h"
#include "BayesianMarkovChain.h"
#include "ContextTreeWeighting.h"
#include "Random.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "DiscreteHiddenMarkovModelOnlineEM.h"
#include "EasyClock.h"
#include "Dirichlet.h"
#include <ctime>
struct ErrorStatistics
{
    std::vector<real> loss;
    std::vector<real> accuracy;
    ErrorStatistics(int T)
    {
        loss.resize(T);
        accuracy.resize(T);
    }
    real total_loss()
    {
        return Sum(loss);
    }
    real total_accuracy()
    {
        return Sum(accuracy);
    }
    void print_result(const char* fname)
    {
        FILE* f = fopen(fname, "w");
        if (!f) return;
        
        int T = loss.size();
        for (int t=0; t<T; ++t) {
            fprintf (f, " %f", loss[t]);
        }
        fprintf (f, "\n");
        fclose (f);
    }
    void print_accuracy(const char* fname)
    {
        FILE* f = fopen(fname, "w");
        if (!f) return;
        
        int T = loss.size();
        for (int t=0; t<T; ++t) {
            fprintf (f, " %f", accuracy[t]);
        }
        fprintf (f, "\n");
        fclose (f);
    }
};




int main (int argc, char** argv)
{
    int n_observations = 4;
    int max_states = 4;
    real prior = 0.5f;
    int T = 100;

    int filename_arg = 1;    
    int max_states_arg = 2;
    int max_horizon_arg = 3;
    int prior_arg = 4;

    setRandomSeed(time(NULL));

    char* fname;
    if (argc < 2) {
        fprintf(stderr, "Usage: bayesian_markov_chain_text filename [max states [max horizon [prior]]]\n");
        exit(-1);
    } else {
        fname = argv[filename_arg];
    }

    if (argc > max_states_arg) {
        max_states = atoi(argv[max_states_arg]);
    }

    if (argc > prior_arg) {
        prior = atof(argv[prior_arg]);
    }


    FILE* file = fopen(fname, "r");
    if (!file) {
        fprintf (stderr, "Error: Could not open file %s\n", fname);
        exit(-1);
    }

    fscanf(file, "%d", &T);

    if (argc > max_horizon_arg) {
        int tmpT = atoi(argv[max_horizon_arg]);
        if (tmpT < T)
            T = tmpT;
    }

    printf("horizon: %d\n", T);
    std::vector<int> data(T);
    n_observations = 0;
    for (int t=0; t<T; ++t) {
        fscanf(file, "%d", &data[t]);
        if (data[t] > n_observations) {
            n_observations = data[t];
        }
        data[t] -= 1;
        //printf("x[%d] = %d\n", t, data[t]);
    }
    fclose(file);

    double bmc_time = 0;
    double bhmc_time = 0;
    double bpsr_time = 0;
    double polya_time = 0;
    double ctw_time = 0;

    double initial_time  = GetCPU();
    double elapsed_time = 0;

    ErrorStatistics bmc_error(T);
    ErrorStatistics bhmc_error(T);
    ErrorStatistics bpsr_error(T);
    ErrorStatistics polya_error(T);
    ErrorStatistics ctw_error(T);


    //logmsg ("Making Bayesian Markov chain\n");
    // our model for the chains
    bool dense = false;
    BayesianMarkovChain bmc(n_observations, 1+max_states, prior, dense);
    BayesianHierarchicalMarkovChain bhmc(n_observations, 1+max_states, prior, dense);
    BVMM bpsr(n_observations, 1+max_states, prior, false, dense);
        //BVMM polya(n_observations, 1+max_states, prior, true, dense);
    DiscreteHiddenMarkovModelOnlineEM polya(1, n_observations);
    ContextTreeWeighting ctw(n_observations, 1+max_states, prior, dense);

    //logmsg ("Making Markov chain\n");
    // the actual chain that generates the data


    bmc.Reset();
    bpsr.Reset();
    ctw.Reset();
    
    
    for (int t=0; t<T; ++t) {
        int observation = data[t];
        
        int bmc_prediction = bmc.predict();
        bmc_error.accuracy[t] += bmc.NextStateProbability(observation);
        if (bmc_prediction != observation) {
            bmc_error.loss[t] += 1.0;
        }

        int bhmc_prediction = bhmc.predict();
        bhmc_error.accuracy[t] += bhmc.NextStateProbability(observation);
        if (bhmc_prediction != observation) {
            bhmc_error.loss[t] += 1.0;
        }
        
        int bpsr_prediction = bpsr.predict();
        bpsr_error.accuracy[t] += bpsr.NextStateProbability(observation);
        if (bpsr_prediction != observation) {
            bpsr_error.loss[t] += 1.0;
        }

        int polya_prediction = polya.predict();
        polya_error.accuracy[t] += (polya.getPrediction())(observation);
//        polya_error.accuracy[t] += polya.NextStateProbability(observation);
        if (polya_prediction != observation) {
            polya_error.loss[t] += 1.0;
        }

        int ctw_prediction = ctw.predict();
        ctw_error.accuracy[t] += ctw.NextStateProbability(observation);
        if (ctw_prediction != observation) {
            ctw_error.loss[t] += 1.0;
        }
        
        //bmc.NextStateProbability(observation);
        //bpsr.NextStateProbability(observation);
        // bmc
        double start_time = GetCPU();
        bmc.ObserveNextState(observation);
        double end_time = GetCPU();
        bmc_time += end_time - start_time;

        // bhmc
        start_time = end_time;
        bhmc.ObserveNextState(observation);
        end_time = GetCPU();
        bhmc_time += end_time - start_time;

        // bpsr
        start_time = end_time;
        bpsr.ObserveNextState(observation);
        end_time = GetCPU();
        bpsr_time += end_time - start_time;

        // polya
        start_time = end_time;
            //polya.ObserveNextState(observation);
        polya.Observe(observation);
        end_time = GetCPU();
        polya_time += end_time - start_time;

        // ctw
        start_time = end_time;
        ctw.ObserveNextState(observation);
        end_time = GetCPU();
        ctw_time += end_time - start_time;
    }
    
    double end_time = GetCPU();
    elapsed_time += end_time - initial_time;
    initial_time = end_time;
        
    
    printf ("# Time -- BMC: %f, BHMC: %f, BPSR: %f, POLYA: %f, CTW: %f\n", 
            bmc_time, bhmc_time, bpsr_time, polya_time, ctw_time);
    printf ("%.1f %.1f %.1f %.1f %.1f # Err\n",
            bmc_error.total_loss(),
            bhmc_error.total_loss(),
            bpsr_error.total_loss(),
            polya_error.total_loss(),
            ctw_error.total_loss());
    printf ("%.1f %.1f %.1f %.1f %.1f # Acc\n",
            bmc_error.total_accuracy(),
            bhmc_error.total_accuracy(),
            bpsr_error.total_accuracy(),
            polya_error.total_accuracy(),
            ctw_error.total_accuracy());

    bmc_error.print_result("bmc_text.error");
    bhmc_error.print_result("bhmc_text.error");
    bpsr_error.print_result("bpsr_text.error");
    polya_error.print_result("polya_text.error");
    ctw_error.print_result("ctw_text.error");

    bmc_error.print_accuracy("bmc_text.accuracy");
    bhmc_error.print_accuracy("bhmc_text.accuracy");
    bpsr_error.print_accuracy("bpsr_text.accuracy");
    polya_error.print_accuracy("polya_text.accuracy");
    ctw_error.print_accuracy("ctw_text.accuracy");
    
}

#endif
