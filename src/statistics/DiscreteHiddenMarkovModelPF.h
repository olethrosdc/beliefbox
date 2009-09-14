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
#ifndef DISCRETE_HIDDEN_MARKOV_MODEL_PF_H
#define DISCRETE_HIDDEN_MARKOV_MODEL_PF_H

#include "DiscreteHiddenMarkovModel.h"
#include "Vector.h"
#include "Dirichlet.h"

#include <vector>

/**
   \ingroup StatisticsGroup
 */
/*@{*/

class DiscreteHiddenMarkovModelPF
{
protected:
    int n_states;
    int n_observations;
    int n_particles;
    std::vector<DiscreteHiddenMarkovModel*> hmm; ///< Transition distribution
    std::vector<DiscreteHiddenMarkovModelStateBelief*> belief; ///< Emission distribution
public:
    Vector P_x;
    Vector log_P_x;
    Vector w;
    Vector log_w;
    std::vector<DirichletDistribution*> state_prior;
    std::vector<DirichletDistribution*> observation_prior;
    DiscreteHiddenMarkovModelPF(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_);
    ~DiscreteHiddenMarkovModelPF();
    real Observe(int x);
    Vector getPrediction();
};


/*@}*/
#endif
