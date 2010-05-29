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
    Matrix _belief; ///< state belief, for EM
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
    void Show();
    real Expectation(std::vector<int>& observations, Matrix& forward_belief, Matrix& backward_belief);
    void Maximisation(Matrix& forward_belief, Matrix& backward_belief);
    real ExpectationMaximisation(std::vector<int>& observations, int n_iterations);
    Matrix& getBelief()
    {
        return _belief;
    }
};


class RandomNumberGenerator;
DiscreteHiddenMarkovModel* MakeRandomDiscreteHMM(int n_states, int n_observations, real stationarity, RandomNumberGenerator* rng);


/** Maintain a belief about the current state
     
    Given a hidden Markov model, maintain a multinomial distribution
    over the states given the history of observations.
    
    More specifically, for a prior \f$\pi\f$ over states and HMM parameters
    \f$\theta\f$, the multinomial distribution \f$B\f$ has parameters
    \f[
    b(i) = \pi(s_t = i | x^t, \theta).
    \f]
      
 */
class DiscreteHiddenMarkovModelStateBelief
{
protected: 
    MultinomialDistribution B; ///< \f$\pi(s_t | x^t)\f$
    int n_states;
public:
    DiscreteHiddenMarkovModel* hmm;
    DiscreteHiddenMarkovModelStateBelief(DiscreteHiddenMarkovModel* hmm_);
    DiscreteHiddenMarkovModelStateBelief(int n_states_);
    DiscreteHiddenMarkovModelStateBelief(int n_states_, int start_state);
    DiscreteHiddenMarkovModelStateBelief(int n_states_, MultinomialDistribution& initial_distribution);
    real Observe(int x);
    real Observe(std::vector<int>& x);
    Vector getBelief();
    Vector getPrediction();
    int predict()
    {
        return ArgMax(getPrediction());
    }
    void Reset();
};


/*@}*/
#endif
