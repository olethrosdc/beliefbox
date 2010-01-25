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
#ifndef BAYESIAN_FMC_H
#define BAYESIAN_FMC_H

#include <vector>
#include <map>
#include "Vector.h"
#include "Matrix.h"
#include "FactoredMarkovChain.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/



/** A mixture of factored Markov Chains.

    @see FactoredMarkovChain
    @see BPSRModel

 */
class BayesianFMC : public FactoredPredictor
{
protected:
	int n_obs;
	int n_actions;
	int n_models;
	real prior;
    std::vector<FactoredMarkovChain*> mc;
    Vector Pr; ///< model probabilities
    Vector logPr; ///< model log probabilities
    Vector Pr_next; ///< observation probabilities
public:
    BayesianFMC (int n_obs_, int n_actions_, int n_models_, real prior_);

    virtual ~BayesianFMC();

    
    /* Training and generation */
    virtual real Observe(int observation);
    virtual real Observe(int action, int observation);
    virtual real ObservationProbability (int action, int observation);
    virtual void Reset();
    virtual int predict(int a);
    
};
/*@}*/
#endif
