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
#include "BayesianMarkovChain.h"
#include "Random.h"
#include "DenseMarkovChain.h"
#include "SparseMarkovChain.h"
#include "EasyClock.h"

struct ErrorStatistics
{
    std::vector<real> loss;
};

void print_result(const char* fname, std::vector<ErrorStatistics>& error)
{
    FILE* f = fopen(fname, "w");
    if (!f) return;
    int max_states = error.size();
    for (int i=0; i<max_states; i++) {
        fprintf (f, "%d", i);
        int T = error[i].loss.size();
        for (int t=0; t<T; ++t) {
            fprintf (f, " %f", error[i].loss[t]);
        }
        fprintf (f, "\n");
    }
    fclose (f);
}


int main (int argc, char** argv)
{
    int n_states = 2;
    int max_states = 8;
    float prior = 1.0f;
    int T = 100;
    int n_iter = 100;

    if (argc>1) {
        T = atoi(argv[1]);
    }

    if (argc>2) {
        n_iter = atoi(argv[2]);
    }

    printf("# Dense posterior \n");
    double dense_time = 0;

    double initial_time  = GetCPU();
    double elapsed_time = 0;

    std::vector<ErrorStatistics> bmc_error(max_states);
    for (int i=0; i<max_states; i++) {
        bmc_error[i].loss.resize(T);
        sparse_bmc_error[i].loss.resize(T);
    }
    for (int iter=0; iter<n_iter; iter++) {
        double remaining_time = (real) (n_iter - iter) * elapsed_time / (real) iter;
        printf ("# iter: %d, %.1f running, %.1f remaining\n", iter, elapsed_time, remaining_time);
        
        for (int i=0; i<max_states; i++) {
            //logmsg ("Making Bayesian Markov chain\n");
            // our model for the chains
            BayesianMarkovChain bmc(n_states, 1+max_states, prior, true);

            //logmsg ("Making Markov chain\n");
            // the actual chain that generates the data
            DenseMarkovChain chain(n_states, i);

            //logmsg ("Creating transitions for Markov chain\n");
            for (int src=0; src<chain.getTotalStates(); ++src) {
                for (int dst=0; dst<n_states; ++dst) {
                    chain.setTransition(src, dst, urandom());
                }
                //chain.setTransition(src, rand()%n_states, 1);
            }

            //logmsg ("Observing chain outputs\n");
            bmc.Reset();
            sparse_bmc.Reset();
            chain.Reset();
        

            for (int t=0; t<T; ++t) {
                int state = chain.generate();

                int bmc_prediction = bmc.generate_static();
                if (bmc_prediction != state) {
                    bmc_error[i].loss[t] += 1.0;
                } else {
                    bmc_error[i].loss[t] += 0.0;
                }

                double start_time = GetCPU();
                bmc.ObserveNextState(state);
                                
                double end_time = GetCPU();
                bmc_time += end_time - start_time;
            }
        }
        double end_time = GetCPU();
        elapsed_time += end_time - initial_time;
        initial_time = end_time;

    }

    printf ("# Time -- Dense: %f, Sparse: %f\n", 
            dense_time, sparse_time);

    print_result("bmc.error", bmc_error);

}

#endif
