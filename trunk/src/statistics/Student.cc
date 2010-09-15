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

#include "Student.h"
#include "SpecialFunctions.h"

/// Default constructor.
///
/// Initialises location to zero and precision to identity.
Student::Student(const int degrees, const int dimension) 
    : n(degrees),
      k(dimension),
      mu(k),
      T(Matrix::Unity(k, k))
{
}

/// Constructor
Student::Student(const int degrees, const Vector& location, const Matrix& precision)
    : n(degrees),
      k(mu.Size()),
      mu(location),
      T(precision)
{
    T.LUDecomposition(det);
}

/// Set the degrees
void Student::setDegrees(const int degrees)
{
    n = degrees;
}
/// Set the location parameter
void Student::setLocation(const Vector& location)
{
    mu = location;
}
/// Set the precision matrix and calculate its determinant
void Student::setPrecision(const Matrix& precision)
{
    T = precision;
    T.LUDecomposition(det);
}

/// Obtain the logarithm of the pdf at x
real Student::log_pdf(const Vector& x) const
{

    real rn = (real) n;
    real rk = (real) k;
    real r = 0.5*(rn + rk);

    Vector delta = x - mu;
    real g = 1 + Mahalanobis2(delta, T, delta) / rn;

    real log_c = logGamma(r) + 0.5 * det - logGamma(0.5*rn) - (0.5*rk)*log(rn*M_PI);
    return log_c - r * log(g);
}

void Student::Show()
{
    printf("n: %d\n", n);
    printf("mu: "); mu.print(stdout);
    printf("T: \n"); T.print(stdout);
        
}
