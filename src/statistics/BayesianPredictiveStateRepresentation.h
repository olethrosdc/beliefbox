/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef BAYESIAN_PSR_H
#define BAYESIAN_PSR_H

#include <vector>
#include <map>
#include "Vector.h"
#include "Matrix.h"
#include "FactoredMarkovChain.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/



/** A tree hierarchy of factored Markov Chains.

    This class uses a context tree to define a probability distribution
    over next observations conditioned on actions.

    @see BVMM
    @see BPSRModel

 */
class BayesianPredictiveStateRepresentation
{
public:
    typedef std::map<int, real, std::greater<int> > BeliefMap;
    typedef BeliefMap::iterator BeliefMapIterator;
protected:
    int total_observations; ///< total number of observations
    int n_obs; ///< number of non-controllable variable states
    int n_actions; ///< number of control variable states
    int n_models; ///< number of models
    Matrix P_obs; ///< Probability of observations for model k
    Matrix Lkoi; ///< Probability of observations for all models up to k
    std::vector<real> weight; ///< temporary weight of model
    std::vector<FactoredMarkovChain*> mc;
    std::vector<real> log_prior;
    Vector Pr; ///< model probabilities
    Vector Pr_next; ///< state probabilities
public:
    std::vector<BeliefMap> beliefs;

    BayesianPredictiveStateRepresentation (int n_obs, int n_actions, int n_models, float prior);

    inline real get_belief_param(int act, int model)
    {
        int src = mc[model]->getContext(act);
        BeliefMapIterator i = beliefs[model].find(src);
		if (i==beliefs[model].end()) {
			return 0.0;
		} else {
			return i->second;
		}
    }

    inline void set_belief_param(int act, int model, real value)
    {
        int src = mc[model]->getContext(act);
        BeliefMapIterator i =  beliefs[model].find(src);
        if (i!=beliefs[model].end()) {
            i->second = value;
        } else {
            beliefs[model].insert(std::make_pair(src, value));
        }
    }

    virtual ~BayesianPredictiveStateRepresentation();

    
    /* Training and generation */
    virtual void Observe(int observation);
    virtual void Observe(int action, int observation);
    virtual real ObservationProbability (int action, int observation);
    virtual void Reset();
    virtual int predict(int a);
    
};
/*@}*/
#endif
