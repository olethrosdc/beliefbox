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


int main (int argc, char** argv)
{
    int n_states = 2;
    int max_states = 10;
    float prior = 1.0f;
    int T = 100;

    if (argc>1) {
        T = atoi(argv[1]);
    }
    for (int i=0; i<max_states; i++) {
        //logmsg ("Making Bayesian Markov chain\n");
        // our model for the chains
        BayesianMarkovChain bmc(n_states, 1+max_states, prior);

        //logmsg ("Making Markov chain\n");
        // the actual chain
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
        chain.Reset();
        for (int t=0; t<T; ++t) {
            int state = chain.generate();
            bmc.ObserveNextState(state);
        }

        printf ("%d", i);
        for (int j=0; j<max_states; ++j) {
            printf (" %f", bmc.Pr[j]);
        }
        printf("\n");
    }
}

#endif
