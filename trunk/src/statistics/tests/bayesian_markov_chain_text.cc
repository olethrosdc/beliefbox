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
#include "BayesianPredictiveStateRepresentation.h"
#include "BayesianMarkovChain.h"
#include "ContextTreeWeighting.h"
#include "Random.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "EasyClock.h"
#include "Dirichlet.h"
#include <ctime>
struct ErrorStatistics
{
    std::vector<real> loss;
    real total_loss()
    {
        return Sum(loss);
    }
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
    float prior = 0.5f;
    int T = 100;
    
    setRandomSeed(time(NULL));

    if (argc > 1) {
        max_states = atoi(argv[1]);
    }


    FILE* file = fopen("data/ifmud.txt", "r");
    fscanf(file, "%d", &T);

    if (argc > 2) {
        int tmpT = atoi(argv[2]);
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
    double bpsr_time = 0;
    double ctw_time = 0;

    double initial_time  = GetCPU();
    double elapsed_time = 0;

    ErrorStatistics bmc_error;
    ErrorStatistics bpsr_error;
    ErrorStatistics ctw_error;
    bmc_error.loss.resize(T);
    bpsr_error.loss.resize(T);
    ctw_error.loss.resize(T);

    //logmsg ("Making Bayesian Markov chain\n");
    // our model for the chains
    bool dense = false;
    BayesianMarkovChain bmc(n_observations, 1+max_states, prior, dense);
    BayesianPredictiveStateRepresentation bpsr(n_observations, 1+max_states, prior, dense);
    ContextTreeWeighting ctw(n_observations, 1+max_states, prior, dense);

    //logmsg ("Making Markov chain\n");
    // the actual chain that generates the data


    bmc.Reset();
    bpsr.Reset();
    //ctw.Reset();
    
    
    for (int t=0; t<T; ++t) {
        int observation = data[t];
        
        int bmc_prediction = bmc.predict();
        if (bmc_prediction != observation) {
            bmc_error.loss[t] += 1.0;
        }
        
        int bpsr_prediction = bpsr.predict();
        if (bpsr_prediction != observation) {
            bpsr_error.loss[t] += 1.0;
        }

        int ctw_prediction = ctw.predict();
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

        // bpsr
        end_time = start_time;
        bpsr.ObserveNextState(observation);
        end_time = GetCPU();
        bpsr_time += end_time - start_time;

        // ctw
        end_time = start_time;
        ctw.ObserveNextState(observation);
        end_time = GetCPU();
        ctw_time += end_time - start_time;
    }
    
    double end_time = GetCPU();
    elapsed_time += end_time - initial_time;
    initial_time = end_time;
        
    
    printf ("# Time -- BHMC: %f, BPSR: %f, CTW: %f\n", 
            bmc_time, bpsr_time, ctw_time);

    printf ("%f %f %f # Error -- BHMC, BPSR, CTW\n", 
            bmc_error.total_loss(),
            bpsr_error.total_loss(),
            ctw_error.total_loss());

    print_result("bmc_text.error", bmc_error);
    print_result("bpsr_text.error", bpsr_error);
    print_result("ctw_text.error", ctw_error);
    
}

#endif
