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

/// This is a generic particle filter for estimating hidden Markov models
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
    Vector getPrediction();
    int predict()
    {
        return ArgMax(getPrediction());
    }
    virtual ~DiscreteHiddenMarkovModelPF();
    virtual real Observe(int x);
    virtual void Reset();
    void Show();
};

/// This particle filter only replaces particles with very small weight
class DiscreteHiddenMarkovModelPF_ReplaceLowest : public DiscreteHiddenMarkovModelPF
{
public:
    real replacement_threshold;
    DiscreteHiddenMarkovModelPF_ReplaceLowest(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_) : 
        DiscreteHiddenMarkovModelPF(threshold, stationarity, n_states_, n_observations_,  n_particles_)
    {
        replacement_threshold = - 2 * log((real) n_particles);
    }  
    virtual ~DiscreteHiddenMarkovModelPF_ReplaceLowest()
    {
    }
    virtual real Observe(int x);
};

/// This particle filter only replaces particles with very small weight
class DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichlet : public DiscreteHiddenMarkovModelPF
{
public:
    real replacement_threshold;
    long T;
    DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichlet(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_) : 
        DiscreteHiddenMarkovModelPF(threshold, stationarity, n_states_, n_observations_,  n_particles_)
    {
        replacement_threshold = - 2 * log((real) n_particles);
        T = 0;
    }  
    virtual ~DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichlet()
    {
    }
    virtual real Observe(int x);
    virtual void Reset()
    {
        T = 0;
    }
};

/// This particle filter only replaces particles with very small weight
class DiscreteHiddenMarkovModelPF_ISReplaceLowest : public DiscreteHiddenMarkovModelPF
{
public:
    real replacement_threshold;
    DiscreteHiddenMarkovModelPF_ISReplaceLowest(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_) : 
        DiscreteHiddenMarkovModelPF(threshold, stationarity, n_states_, n_observations_,  n_particles_)
    {
        replacement_threshold = - 2 * log((real) n_particles);
    }  
    virtual ~DiscreteHiddenMarkovModelPF_ISReplaceLowest()
    {
    }
    virtual real Observe(int x);
};

/// This particle filter only replaces some particles but uses all history to 
class DiscreteHiddenMarkovModelPF_ReplaceLowestExact: public DiscreteHiddenMarkovModelPF
{
public:
    real replacement_threshold;
    std::vector<int> history; ///< observation history
    DiscreteHiddenMarkovModelPF_ReplaceLowestExact(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_) : 
        DiscreteHiddenMarkovModelPF(threshold, stationarity, n_states_, n_observations_,  n_particles_)
    {
        replacement_threshold = - 2 * log((real) n_particles);
    }  
    virtual ~DiscreteHiddenMarkovModelPF_ReplaceLowestExact()
    {
    }
    virtual real Observe(int x);
};


/// This particle replaces all particles 
class DiscreteHiddenMarkovModelRBPF : public DiscreteHiddenMarkovModelPF
{
protected:
    long t;
public:
    DiscreteHiddenMarkovModelRBPF(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_) : 
        DiscreteHiddenMarkovModelPF(threshold, stationarity, n_states_, n_observations_,  n_particles_), t(0)
    {
    }  
    virtual ~DiscreteHiddenMarkovModelRBPF()
    {
    }
    virtual real Observe(int x);
    virtual void Reset();
};




/** This is a mixture of discrete HMM particle filters */
template <class DHMM_Filter>
class DHMM_PF_Mixture
{
protected:
    int n_observations;
    int n_particles;
    int max_states;
    FiniteMixture<DHMM_Filter, int, Vector> mixture;
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
            DHMM_Filter* hmm = new DHMM_Filter(threshold, stationarity, i + 1, n_observations, n_particles);
            mixture.SetComponent(hmm, i);
        }
    }
    
    ~DHMM_PF_Mixture()
    {
        for (int i=0; i<max_states; ++i) {
            delete mixture.GetComponent(i);
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
    void Reset()
    {
        for (int k=0; k<max_states; ++k) {
            mixture.GetComponent(k)->Reset();
        }
    }
};

/*@}*/
#endif
