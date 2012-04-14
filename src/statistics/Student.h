/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef STUDENT_DISTRIBUTION_H
#define STUDENT_DISTRIBUTION_H

#include "Distribution.h"
#include "Matrix.h"
#include "Vector.h"

class MultivariateNormal;

/**
   \ingroup StatisticsGroup
 */
/*@{*/

/** \brief Multivariate Student t-distribution.
    
    A multivariate t distribution with \f$n\f$ degrees of freedom,
    location vector \f$\mu \in R^k\f$ and precision matrix
    \f$T \in R^{k \times k}\f$ has a density defined for all \f$x \in R^k\f$.
    \f[
    f(x \mid n, \mu, T) = c\left[1 + \frac{1}{n} (x-\mu)'T(x-\mu)\right]^{-(n+k)/2},
    \f]
    where
    \f[
    c = \frac{\Gamma[(n + k) / 2] |T|^{1/2}}{\Gamma(n / 2)(n \pi)^{k/2}}.
    \f]
*/
class Student
{
private:
    MultivariateNormal* sampler;
public:
    int n; ///< Degrees of freedom
    const int k; ///< Dimensionality
    Vector mu; ///< location
    Matrix T; ///< precision
    real det; ///< precision determinant
    Student(const int dimension);
    Student(const int degrees, const Vector& location, const Matrix& precision);
    virtual ~Student();
    void setDegrees(const int degrees);
    void setLocation(const Vector& location);
    void setPrecision(const Matrix& precision);
    real log_pdf(const Vector& x) const;
    /// Obtain the pdf at x
    real pdf(const Vector& x) const
    {
        return exp(log_pdf(x));
    }
    void Show() const;
    Vector generate() const;
};

/*@}*/

#endif
