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
#ifndef BVMM_H
#define BVMM_H

#include "BayesianMarkovChain.h"
#include <vector>
#include <map>
#include "Vector.h"
#include "Matrix.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/



/// A Markov Chain
class BVMM : public BayesianMarkovChain
{
protected:
    bool polya;
    Matrix P_obs; ///< Probability of observations for model k
    Matrix Lkoi; ///< Probability of observations for all models up to k
    std::vector<real> weight; ///< temporary weight of model
    typedef std::map<int, real, std::greater<int> > BeliefMap;
    typedef BeliefMap::iterator BeliefMapIterator;
public:
    std::vector<BeliefMap> beliefs;

    BVMM (int n_states, int n_models, float prior, bool polya_, bool dense = false);

    inline real get_belief_param(int model)
    {
        int src = mc[model]->getCurrentState();
        BeliefMapIterator i = beliefs[model].find(src);
		if (i==beliefs[model].end()) {
			return 0.0;
		} else {
			return i->second;
		}
    }

    inline void set_belief_param(int model, real value)
    {
        int src = mc[model]->getCurrentState();
        BeliefMapIterator i =  beliefs[model].find(src);
        if (i!=beliefs[model].end()) {
            i->second = value;
        } else {
            beliefs[model].insert(std::make_pair(src, value));
        }
    }

    virtual ~BVMM();

    
    /* Training and generation */
    virtual void ObserveNextState (int state);
    inline void Observe(int x) {
        ObserveNextState(x);
    }
    virtual real NextStateProbability (int state);
    virtual void Reset();
    virtual int generate();
    virtual int predict();
    
};
/*@}*/
#endif
