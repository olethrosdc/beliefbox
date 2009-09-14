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
#ifndef DISCRETE_HIDDEN_MARKOV_MODEL_H
#define DISCRETE_HIDDEN_MARKOV_MODEL_H

#include "DiscreteHiddenMarkovModel.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/

class DiscreteHiddenMarkovModelPF
{
protected:
    int n_particles;
    std::vector<DiscreteHiddenMarkovModel*> hmm; ///< Transition distribution
    std::vector<DiscreteHiddenMarkovModelStateBelief*> belief; ///< Emission distribution
    std::vector<real> w;
    std::vector<real> log_w;
public:
    
};


/*@}*/
#endif
