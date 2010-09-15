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

/// Multivariate student t-distribution
class Student
{
 protected:
    int n; ///< Degrees of freedom
    const int k; ///< Dimensionality
    Vector mu; ///< location
    Matrix T; ///< precision
    real det; ///< precision determinant
 public:
    Student(const int degrees, const int dimension);
    Student(const int degrees, const Vector& location, const Matrix& precision);
    void setDegrees(const int degrees);
    void setLocation(const Vector& location);
    void setPrecision(const Matrix& precision);
    real log_pdf(const Vector& x) const;
    /// Obtain the pdf at x
    real pdf(const Vector& x) const
    {
        return exp(log_pdf(x));
    }
    void Show();
};

#endif
