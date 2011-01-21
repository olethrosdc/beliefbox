/* -*- Mode: C++ -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MOMENT_MATCHING_BETA_ESTIMATE_H
#define MOMENT_MATCHING_BETA_ESTIMATE_H

#include "BetaDistribution.h"

/** A multivariate beta-product distribution with moment matching
     */
class MomentMatchingBetaEstimate
{
protected:
    std::vector<BetaDistribution*> beta;
    int n_dim;
public:
    int T;
    Vector a; ///< lower bound
    Vector b; ///< upper bound
    Vector d; ///< difference
    Vector c; ///< inverse difference
    MomentMatchingBetaEstimate(const Vector& lower_bound, const Vector& upper_bound);
    void Reset();
    virtual ~MomentMatchingBetaEstimate();
    virtual Vector generate();
    virtual Vector generate() const;
    /// Note that this the marginal likelihood!
    virtual real pdf(const Vector& x) const;
    /// The marginal log-likelihood
    virtual real log_pdf(const Vector& x) const;
    virtual const Vector& getMean() const;
    real Observe(const Vector& x);
    void Show();
    Vector transform(const Vector& x) const;
    Vector inverse_transform(const Vector& x) const;
};


#endif
