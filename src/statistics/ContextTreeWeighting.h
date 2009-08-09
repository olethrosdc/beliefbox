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
#ifndef CONTEXT_TREE_WEIGHTING_H
#define CONTEXT_TREE_WEIGHTING_H

#include "BayesianMarkovChain.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/


/// A Markov Chain
class ContextTreeWeighting : public BayesianMarkovChain
{
public:
    int n_observations;

    ContextTreeWeighting (int n_states, int n_models, float prior, bool dense = false);


    virtual ~ContextTreeWeighting();

    
    /* Training and generation */
    virtual void ObserveNextState (int state);
    virtual float NextStateProbability (int state);
    virtual void Reset();
    virtual int generate();
    virtual int predict();
    
};
/*@}*/
#endif
