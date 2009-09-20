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
#ifndef FINITE_MIXTURE_H
#define FINITE_MIXTURE_H

#include <vector>
#include "Vector.h"
#include "real.h"

/**
   \ingroup StatisticsGroup
*/
/*@{*/

/** A finite mixture template
    
   T is the type of the distribution
   X is the data type
   W is the type of distribution

   T is required to implement two functions:
   
   real Observe(X x);
   
   The argument of this function is the currently
   observed data, while the return value is the likelihood
   of the data under the model implemented by T. 
   
   After (and only after) calculating the likelihood, the model should be
   adapted to take into account the new data.

   W getPrediction()

   This should return a distribution on the current data.
*/
template <typename T, typename X, typename W>
class FiniteMixture
{
protected:
    int n_components; ///< number of mixture components
    std::vector<T*> Q; ///< pointers to mixture components
    Vector log_Pr; ///< log probabiility of components
    Vector Pr; ///< probability of components
    Vector log_L; ///< log likelihood
public:

/// Construct an empty model
    FiniteMixture(int n) : n_components(n),
                           Q(n_components),
                           log_Pr(n_components),
                           Pr(n_components),
                           log_L(n_components)
    {
        assert (n_components > 0);
        ResetPriorToUniform();
    }
    
    /// Construct a model from a vector of components
    FiniteMixture(std::vector<T*>& Q_) : n_components(Q_.size()),
                                         Q(Q_),
                                         log_Pr(n_components),
                                         Pr(n_components),
                                         log_L(n_components)
    {
        assert (n_components > 0);
        ResetPriorToUniform();
    }

    ~FiniteMixture()
    {
#if 0
        for (int k=0; k<n_components; ++k) {
            printf ("%f ", Pr[k]);
        }
        printf("# final weight\n");
#endif
    }
    
    /// Reset to a uniform prior
    void ResetPriorToUniform()
    {
        real prior = 1.0 / (real) n_components;
        real log_prior = log(prior);
        for (int k=0; k<n_components; k++) {
            Pr[k] = prior;
            log_Pr[k] = log_prior;
        }
    }

    /// Make a specific prior
    void SetPrior(Vector P)
    {
        for (int k=0; k<n_components; k++) {
            Pr[k] = P[k];
            log_Pr[k] = log(Pr[k]);
        }
    }

    /// Set the model for a component
    void SetComponent(T* q, int k)
    {
        assert(k >= 0 && k < n_components);
        assert(q);
        Q[k] = q;
    }

    /// Get the model of a component
    T* GetComponent(int k)
    {
        assert(k >= 0 && k < n_components);
        return Q[k];
    }
    
    /// Observe new data, report likelihood and adapt prior over models.
    real Observe(X x)
    {
        real log_sum = LOG_ZERO;

        // calculate log p(x,k) = log p(x|k)+ log p(k) for each k
        // and log p(x) = log sum_k p(x,k)
        for (int k=0; k<n_components; ++k) {
            log_L[k] = log(Q[k]->Observe(x)) + log_Pr[k];
            log_sum = logAdd(log_sum, log_L[k]);
        }

        // log p(k|x) = log p(x,k) - log  p(x)
        log_Pr = log_L - log_sum;
        Pr = exp(log_Pr);
        assert(approx_eq(Pr.Sum(), 1));
        return exp(log_sum);
    }

    /// Get the prediction of the k-th component
    W getPrediction(int k)
    {
        return Q[k]->getPrediction();
    }

    /// Get the weight of the
    real getWeight(int k)
    {
        return Pr[k];
    }

    
    
};
/*@}*/

#endif
