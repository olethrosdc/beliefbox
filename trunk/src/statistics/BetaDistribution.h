/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BETA_DISTRIBUTION_H
#define BETA_DISTRIBUTION_H

#include "Distribution.h"

class BetaDistribution : public ConjugatePrior
{
public:
    real alpha; ///< number of ones obsered
    real beta; ///< number of zeros observed

    /// Beta distribution default constructor: alpha,beta=1, the
    /// uniform distribution
    BetaDistribution()
    {
        alpha = 1.0f;
        beta = 1.0f;
    }
    /// Give your own alpha/beta values
    BetaDistribution(real alpha, real beta)
    {
        this->alpha = alpha;
        this->beta = beta;
    }
    BetaDistribution(const BetaDistribution& rhs)
    {
        alpha = rhs.alpha;
        beta = rhs.beta;
    }
    BetaDistribution& operator = (const BetaDistribution & rhs)
    {
        if (this != &rhs) {
            alpha = rhs.alpha;
            beta = rhs.beta;
        }
        return *this;
    }

    virtual real marginal_pdf(real x) const;
    virtual real pdf(real x) const;
    real log_pdf(real x) const;
    virtual void calculatePosterior(real x);
    virtual void setMean(real mean);
    virtual void setVariance(real var);
    virtual real getMean() const;
    virtual real getVariance(); 
    virtual real generate();
	real generate() const;
    virtual real generateMarginal();
    real Observe(real x);
    real setMaximumLikelihoodParameters(const std::vector<real>& x,
                                        int n_iterations);
};



/// Conjugate prior to gamma distribution with unknown shape and
/// inverse scale
class BetaDistributionMCPrior
{
public:
    real kappa; ///< exponential parameter 1
    int lambda; ///< exponential parameter 2
    BetaDistributionMCPrior(real kappa_, real lambda_); ///< constructor
    virtual ~BetaDistributionMCPrior(); ///< destructor
    virtual real LogLikelihood(std::vector<real>& x, int K) const; 
};


#endif
