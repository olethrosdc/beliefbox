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

// a finite mixture template
//
// T is the type of the distribution
// X is the data type
template <typename T, typename X>
class FiniteMixture
{
protected:
    int n_components; ///< number of mixture components
    std::vector<T*> Q; ///< pointers to mixture components
    Vector log_Pr; ///< log probabiility of components
    Vector Pr ///< probability of components
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
        ResetPriorToUniform()
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
        assert(q);
        Q[k] = q;
    }
    
    real Observe(X x)
    {
        real log_sum = LOG_ZERO;

            // calculate log p(x,k) = log p(x|k)+ log p(k) for each k
            // and log p(x) = log sum_k p(x,k)
        for (int k=0; k<n_components; ++k) {
            log_L[k] = log(Q[k]->Observe(x)) + log_P[k];
            log_sum = logAdd(log_sum, log_L[k]);
       }

            // llog p(k|x) = log p(x,k) - log  p(x)
        log_P = log_L - log_sum;
        P = exp(log_P);
        
    }
    
};
/*@}*/

#endif
