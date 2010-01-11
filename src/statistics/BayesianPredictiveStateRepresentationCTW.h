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
#ifndef BAYESIAN_PSR_CTW_H
#define BAYESIAN_PSR_CTW_H

#include "BayesianPredictiveStateRepresentation.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/



/** A tree hierarchy of factored Markov Chains.

    This class uses a context tree to define a probability distribution
    over next observations conditioned on actions.

    @see BVMM
    @see BPSRModel
    @see BayesianPredictiveStateRepresentation
 */
class BayesianPredictiveStateRepresentationCTW : public BayesianPredictiveStateRepresentation
{
public:
    BayesianPredictiveStateRepresentationCTW (int n_obs, int n_actions, int n_models, float prior);

    virtual ~BayesianPredictiveStateRepresentationCTW();

    
    /* Training and generation */
    virtual void Observe(int observation);
    virtual void Observe(int action, int observation);
    virtual real ObservationProbability (int action, int observation);
    virtual int predict(int a);

};
/*@}*/
#endif
