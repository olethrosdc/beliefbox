/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef WISHART_H
#define WISHART_H

#include "Distribution.h"

/// Wishart probability distribution
class Wishart : public AbstractDistribution<Matrix>
{
public:
    int k; ///< dimensionality
    real n; ///< degrees of freedom
    Matrix V; ///< scale matrix
    Wishart();
    Wishart(real& n_, Matrix& V_);
    virtual ~Wishart();
    virtual void generate(Matrix& x) const;
    virtual Matrix generate() const;
    virtual real pdf(const Matrix& x) const =0;
    virtual real log_pdf(const Matrix& x) const
    {
        return log(pdf(x));
    }

};

#endif
