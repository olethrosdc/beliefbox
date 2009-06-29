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
#include "Random.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "EasyClock.h"
#include "Dirichlet.h"
#include <ctime>
struct ErrorStatistics
{
    std::vector<real> loss;
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
    real stationarity = 0.99;
    
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

    double bmc_time = 0;
    double bpsr_time = 0;

    double initial_time  = GetCPU();
    double elapsed_time = 0;

    ErrorStatistics oracle_error;
    ErrorStatistics bmc_error;
    ErrorStatistics bpsr_error;
    oracle_error.loss.resize(T);
    bmc_error.loss.resize(T);
    bpsr_error.loss.resize(T);

    for (int iter=0; iter<n_iter; iter++) {
        double remaining_time = (real) (n_iter - iter) * elapsed_time / (real) iter;
        printf ("# iter: %d, %.1f running, %.1f remaining\n", iter, elapsed_time, remaining_time);
        
        
        //logmsg ("Making Bayesian Markov chain\n");
        // our model for the chains
        bool dense = false;
        BayesianMarkovChain bmc(n_observations, 1+max_states, prior, dense);
        
        BayesianPredictiveStateRepresentation bpsr(n_observations, 1+max_states, prior, dense);

        //logmsg ("Making Markov chain\n");
        // the actual chain that generates the data
        DenseMarkovChain chain(n_mc_states, 1);
        DirichletDistribution dirichlet(n_mc_states, 1.0);
        std::vector<std::vector<real> >X(n_mc_states);
        for (int i=0; i<n_mc_states; ++i) {
            X[i].resize(n_observations);
            real sum = 0.0;
            for (int j=0; j<n_observations; ++j) {
                X[i][j] = 0.1*true_random(false);
                if (i==j) {
                    X[i][j] += 1.0;
                }
                sum += X[i][j];
            }
            for (int j=0; j<n_observations; ++j) {
                X[i][j] /=  sum;
            }

        }
        //logmsg ("Creating transitions for Markov chain\n");
        for (int src=0; src<n_mc_states; ++src) {
            Vector Pr(n_mc_states);
            for (int i=0; i<n_mc_states; ++i) {
                Pr[i] = exp(10.0*true_random(false));
            }
            Pr /= Pr.Sum();
            Pr *= (1 - stationarity);
            Pr[src] += stationarity;
            Pr /= Pr.Sum();
            real sum = 0.0;
            for (int dst=0; dst<n_observations; ++dst) {
                //printf ("P(%d|%d)=%f\n", dst, src, Pr[dst]);
                chain.setTransition(src, dst, Pr[dst]);
                sum += Pr[dst];
            }
            //printf("  Tot: %f\n", sum);
        }

        //logmsg ("Observing chain outputs\n");
        bmc.Reset();
        bpsr.Reset();
        chain.Reset();
        

        for (int t=0; t<T; ++t) {
            int predicted_state = chain.GenerateStatic();
            int predicted_observation = DiscreteDistribution::generate(X[predicted_state]);
            int state = chain.generate();
            int observation = DiscreteDistribution::generate(X[state]);
            //printf("%d #obs\n", observation);

            //printf("%d |BHMC| ", state);

            if (predicted_observation != observation) {
                oracle_error.loss[t] += 1.0;
            }

            int bmc_prediction = bmc.predict();
            if (bmc_prediction != observation) {
                bmc_error.loss[t] += 1.0;
            }

            //printf("  |BPSR| ");
            int bpsr_prediction = bpsr.predict();
            if (bpsr_prediction != observation) {
                bpsr_error.loss[t] += 1.0;
            }
                
            //bmc.NextStateProbability(observation);
            //bpsr.NextStateProbability(observation);
            double start_time = GetCPU();
            bmc.ObserveNextState(observation);
                                
            double end_time = GetCPU();
            bmc_time += end_time - start_time;

            end_time = start_time;
            bpsr.ObserveNextState(observation);
                                
            end_time = GetCPU();
            bpsr_time += end_time - start_time;
        }

        double end_time = GetCPU();
        elapsed_time += end_time - initial_time;
        initial_time = end_time;
        
    }

    printf ("# Time -- BHMC: %f, BPSR: %f\n", 
            bmc_time, bpsr_time);

    real inv_iter = 1.0 / (real) n_iter;
    for (int t=0; t<T; ++t) {
        oracle_error.loss[t] *= inv_iter;
        bmc_error.loss[t] *= inv_iter;
        bpsr_error.loss[t] *= inv_iter;
    }
    print_result("oracle.error", oracle_error);
    print_result("bmc.error", bmc_error);
    print_result("bpsr.error", bpsr_error);

}

#endif
