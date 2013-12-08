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
  //det = T.det();
  T.LUDecomposition(det);
}

/// Constructor
Student::Student(const int degrees, const Vector& location, const Matrix& precision)
    : sampler(new MultivariateNormal(location.Size())),
      n(degrees),
      k(location.Size()),
      mu(location),
      T(precision)
{
  //det = T.det();
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
	//det = T.det();
    //printf("New Precision det:%f\n", det);
    //T.print(stdout);
}

/** Obtain the logarithm of the pdf at \f$x \in R^k\f$
	
	\f[
	f(x, T, d) = 
	\frac{\Gamma((d+k)/2)}{\Gamma(d/2)}
	|T|^{1/2} (d\pi)^{k/2} 
	\left(1 + x^\top T x / d
	\right)^{-(d+k)/2}
	\f]
 */
real Student::log_pdf(const Vector& x) const
{
	Vector delta = x - mu;
    real degree = (real) n;
    real g = 1.0 + Mahalanobis2(delta, T, delta) / degree;
//	printf("degree = %f k = %f  r = %f\n",degree, (real)(k), 0.5*(degree + k));
    real l1 = logGamma(0.5 * ( degree + (real) k));
	real l2 = - logGamma(0.5 * degree);
	real l3 = + 0.5 * log(det);
	real l4 = - (0.5 * (real) k)*log(degree * M_PI); 
	real l5 =  - 0.5 * (degree + k) * log(g);
    real log_p = l1 + l2 + l3 + l4 + l5;
//	printf ("g: %f l1: %f l2: %f l3: %f l4: %f l5: %f log_p: %f\n", 
//		g, l1, l2, l3, l4, l5, log_p);
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
