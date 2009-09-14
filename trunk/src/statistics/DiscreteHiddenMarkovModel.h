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

#include "HiddenMarkovModel.h"
#include "MultinomialDistribution.h"
#include "Matrix.h"
#include <vector>

/**
   \ingroup StatisticsGroup
 */
/*@{*/

class DiscreteHiddenMarkovModel
{
protected:
    int n_states;
    int n_observations;
    std::vector<MultinomialDistribution> P_S; ///< Transition distribution
    std::vector<MultinomialDistribution> P_X; ///< Emission distribution
    int current_state;
public:
    DiscreteHiddenMarkovModel(Matrix& Pr_S, Matrix& Pr_X);
    DiscreteHiddenMarkovModel(int n_states_, int n_observations_);
    void Reset();
    virtual ~DiscreteHiddenMarkovModel ();
    int getCurrentState()
    {
        return current_state;
    }
    virtual int generate();
    virtual int generate_static();
    int getNStates()
    {
        return n_states;
    }
    int getNObservations()
    {
        return n_observations;
    }
    std::vector<MultinomialDistribution>& getStateProbablities()
    {
        return P_S;
    }
    std::vector<MultinomialDistribution>& getObservationProbablities()
    {
        return P_X;
    }
        /// Given a source state, get the probability of a destination state
    inline real& PrS(int src, int dst)
    {
        return P_S[src].Pr(dst);
    }
        /// Get probability of observation given source state
    inline real& PrX(int s, int x)
    {
        return P_X[s].Pr(x);
    }
};





class DiscreteHiddenMarkovModelStateBelief
{
protected:
    MultinomialDistribution B;
    int n_states;
public:
    DiscreteHiddenMarkovModel* hmm;
    DiscreteHiddenMarkovModelStateBelief(DiscreteHiddenMarkovModel* hmm_);
    DiscreteHiddenMarkovModelStateBelief(int n_states_);
    DiscreteHiddenMarkovModelStateBelief(int n_states_, int start_state);
    DiscreteHiddenMarkovModelStateBelief(int n_states_, MultinomialDistribution& initial_distribution);
    real Observe(int x);
    Vector getBelief();
    Vector getPrediction();
};


/*@}*/
#endif
