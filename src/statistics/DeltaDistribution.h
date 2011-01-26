/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#include "real.h"
#include "Vector.h"

/** A delta distribution with a uniform prior
 */
class DeltaUniformDistribution
{
public: 
    DeltaUniformDistribution(const Vector& lower_bound_,
                             const Vector& upper_bound_);
    real pdf(const Vector& x) const;
    real Observe(const Vector& x);
protected:
    int n_dim; ///< the number of dimensions
    Vector lower_bound; ///< lower bound on where it may lie
    Vector upper_bound; ///< upper bound on where it may lie
    Vector mean; ///< mean
    int T; ///< number of observations

    
};
