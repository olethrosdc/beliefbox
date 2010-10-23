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
#include "Distribution.h"
#include <vector>

/** Gaussian process. Not implemented!
    
    TODO Implement gaussian process
 */
class GaussianProcess : public VectorDistribution
{
protected:
    std::vector<Vector> kernels; ///< list of kernels
    int N; ///< number of kernels
    int d; ///< number of dimensions
public:
    GaussianProcess(int d_);
    GaussianProcess(std::vector<Vector> k);
    virtual ~GaussianProcess();
    virtual Vector generate();
    virtual real pdf(Vector x);
    virtual void update(Vector* x);
};

#endif

