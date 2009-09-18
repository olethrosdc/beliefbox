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
#include "FiniteMixture.h"
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
    int predict()
    {
        return ArgMax(getPrediction());
    }
    void Reset();
};

/** This is a mixture of discrete HMM particle filters */
class DHMM_PF_Mixture
{
protected:
    int n_observations;
    int n_particles;
    int max_states;
    FiniteMixture<DiscreteHiddenMarkovModelPF, int, Vector> mixture;
public:
    DHMM_PF_Mixture(real threshold,
                    real stationarity,
                    int n_observations_,
                    int n_particles_,
                    int max_states_) : n_observations(n_observations_),
                                       n_particles(n_particles_),
                                       max_states(max_states_),
                                       mixture(max_states)
    {
        Vector P(max_states);
        for (int i=0; i<max_states; ++i) {
            P[i] = exp(- (real) i * (i + n_observations));
        }
        P /= P.Sum();
        mixture.SetPrior(P);
        for (int i=0; i<max_states; ++i) {
            DiscreteHiddenMarkovModelPF* hmm = new DiscreteHiddenMarkovModelPF(threshold, stationarity, i + 1, n_observations, n_particles);
            mixture.SetComponent(hmm, i);
        }
    }
    
    real Observe(int x)
    {
        return mixture.Observe(x);
    }

    Vector getPrediction()
    {
        Vector p(n_observations);
        for (int k=0; k<max_states; ++k) {
            p += mixture.getPrediction(k) * mixture.getWeight(k);
        }
        return p;
    }
    
    Vector getWeights()
    {
        Vector p(max_states);
        for (int k=0; k<max_states; ++k) {
            p[k] = mixture.getWeight(k);
        }
        return p;
    }

    int predict()
    {
        return ArgMax(getPrediction());
    }
};

/*@}*/
#endif
