/* -*- Mode: c++;  -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef SPARSE_MARKOVCHAIN_H
#define SPARSE_MARKOVCHAIN_H

#include <vector>
#include "real.h"
#include "MarkovChain.h"
#include "Vector.h"

#include "SparseTransitions.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/

/** A sparse implementation of a Markov chain
 */
class SparseMarkovChain : public MarkovChain
{
protected:
    int n_transitions; ///< total number of transitions
	SparseTransitions transitions; ///< history-wide transition table

    real threshold;
public:
    SparseMarkovChain (int n_states, int mem_size);
    virtual ~SparseMarkovChain ();

    /* probabilities */
    virtual real getTransition (MCState src, int dst);
    virtual real getProbability (MCState src, int dst);
    virtual void getProbabilities(MCState src, std::vector<real>& p);
    virtual void getNextStateProbabilities(std::vector<real>& p);
    virtual real pdf(MCState src, Vector q);
    virtual void setTransition (MCState src, int dst, real value);
    virtual void setThreshold (real threshold);


    /* Training and generation */
    virtual int PushState (int state);
    virtual real ObserveNextState (int state);
    virtual real NextStateProbability (int state);
    virtual void Reset();
    virtual int GenerateStatic ();
    virtual int generate ();
    
    /* Debug functions */
    virtual int ShowTransitions ();
};
/*@}*/
#endif
