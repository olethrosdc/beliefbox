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
#include "MultivariateNormal.h"
#include "ranlib.h"

/// Default constructor.
///
/// Initialises location to zero and precision to identity.
Student::Student(const int dimension) 
    : sampler(new MultivariateNormal(dimension)),
      n(1),
      k(dimension),
      mu(k),
      T(Matrix::Unity(k, k))
{
}

/// Constructor
Student::Student(const int degrees, const Vector& location, const Matrix& precision)
    : sampler(new MultivariateNormal(mu.Size())),
      n(degrees),
      k(mu.Size()),
      mu(location),
      T(precision)
{
    T.LUDecomposition(det);
}

Student::~Student()
{
    delete sampler;
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
    //printf("New Precision det:%f\n", det);
    //T.print(stdout);
}

/** Obtain the logarithm of the pdf at x.   
 */
real Student::log_pdf(const Vector& x) const
{
    Vector delta = x - mu;
    real degree = (real) n;
    real g = 1 + Mahalanobis2(delta, T, delta) / degree;
    real log_c = logGamma(0.5 * (degree + k))
        + 0.5 * log(det)
        - logGamma(0.5 * degree)
        - (0.5 * k)*log(n*M_PI); // something is maybe missing here
    real log_p = log_c - (0.5 * (degree + k)) * log(g);
    return log_p;
}

void Student::Show() const
{
    printf("# Student parameters\n");
    printf("#   degrees: %d\n", n);
    printf("#  location: "); mu.print(stdout);
    printf("# precision: "); T.print(stdout);
    printf("# --\n");
        
}

/** Generate a sample.

    Simply draw use a normal and a chi^2 variate dude!
*/
Vector Student::generate() const
{
    sampler->setAccuracy(T);
    Vector v = sampler->generate();
    real z = genchi((real) n);
    //v.print(stdout);
    //printf ("%f -> ", z);
    v /= z;
    //v.print(stdout);
    return mu + v;
}
