/* -*- Mode: c++;  -*- */
/*VER: $Id: MarkovChain.h,v 1.7 2006/08/17 05:35:00 olethros Exp $*/
// copyright (c) 2004 by Christos Dimitrakakis <dimitrak@idiap.ch>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef BAYESIAN_MARKOVCHAIN_H
#define BAYESIAN_MARKOVCHAIN_H

#include "MarkovChain.h"
#include <vector>
#include "Vector.h"

/**
   \ingroup StatisticsGroup
*/
/*@{*/

/// A Markov Chain
class BayesianMarkovChain
{
public:
    int n_states; ///< number of distinct states
    int n_models; ///< number of models
    std::vector<MarkovChain*> mc; ///< Markov chain
    std::vector<real> log_prior;
    Vector Pr; ///< model probabilities
    Vector Pr_next; ///< state probabilities
    int n_observations;

    BayesianMarkovChain (int n_states, int n_models, real prior, bool dense = false);
    virtual ~BayesianMarkovChain();


    /* Training and generation */
    virtual void ObserveNextState (int state);
    virtual real NextStateProbability (int state);
    virtual void Reset();
    virtual int generate();
    virtual int predict();
    
    /* Helper functions */  
    int CalculateStateID ();
    int  PushState (int state);
    
    /* Debug functions */
    void ShowTransitions ();
};
/*@}*/
#endif
