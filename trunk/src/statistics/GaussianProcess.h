/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GAUSSIAN_PROCESS_H
#define GAUSSIAN_PROCESS_H

#include "Vector.h"
#include "Matrix.h"
#include "Distribution.h"
#include <vector>

/** Gaussian process. 
    
    This is a {\em conditional} distribution.
 */
class GaussianProcess
{
protected:
    Matrix Sigma_p;
    Matrix Accuracy;
    Matrix A;
    real noise_variance;
    Matrix X2; ///< observation co-variance.
    Vector mean;
    Matrix covariance;
public:
    GaussianProcess(Matrix& Sigma_p_,
                    real noise_variance_);
    virtual ~GaussianProcess();
    virtual Vector generate();
    virtual real pdf(Vector& x, real y);
    virtual void Observe(Vector& x, real y);
};

#endif

