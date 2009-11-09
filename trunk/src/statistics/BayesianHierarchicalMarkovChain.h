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
#ifndef BAYESIAN_HIERARCHICAL_MARKOVCHAIN_H
#define BAYESIAN_HIERARCHICAL_MARKOVCHAIN_H

#include "BayesianMarkovChain.h"
#include "Matrix.h"

/**
   \ingroup StatisticsGroup
*/
/*@{*/

/// A hierarchical Markov Chain
class BayesianHierarchicalMarkovChain : public BayesianMarkovChain
{
public:
    Vector weight; ///< weights
    Vector log_weight; ///< log weights
    Matrix Lkoi; ///< P(x|B_k) = P(x|B_{k-1}) (1 - w_k) + w_k P(x|M_k)
    Matrix P_obs; ///< P(x|M_k) 
    Vector P_next;
    int prediction;
    BayesianHierarchicalMarkovChain (int n_states, int n_models, float prior, bool dense = false);
    virtual ~BayesianHierarchicalMarkovChain();

    /* Training */
    virtual void ObserveNextState (int state);
    virtual int predict();
    virtual real NextStateProbability(int state)
    {
        return Pr_next[state];
    }
};
/*@}*/
#endif
