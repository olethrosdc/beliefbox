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

#include "Dirichlet.h"
#include "ranlib.h"
#include "SpecialFunctions.h"

/// Create a placeholder Dirichlet
DirichletDistribution::DirichletDistribution()
	: n_observations(0)
{
    n = 0;
    alpha_sum = 1.0;
}

/// Create a Dirichlet with uniform parameters
DirichletDistribution::DirichletDistribution(int n, real p)
	: n_observations(0)
{
    resize(n, p);
	alpha_sum = n * p;
}

/// Initialise parameters from a vector
DirichletDistribution::DirichletDistribution(const Vector& x) : n(x.Size()), alpha(x), n_observations(0)
{
    for (int i=0; i<n; ++i) {
        assert(alpha(i) >= 0);
    }
    alpha_sum = alpha.Sum();
}

/// Resize the dirichlet distribution, setting all parameters to p.
void DirichletDistribution::resize(int n, real p)
{
    this->n = n;
    alpha.Resize(n);
    for (int i=0; i<n; ++i) {
        alpha(i) = p;
    }
    alpha_sum = (real) n * p;
}


/// Destructor
DirichletDistribution::~DirichletDistribution()
{
}

/// Generate a multinomial vector
Vector DirichletDistribution::generate() const
{
    Vector x(n);
    generate(x);
    return x;
}

/// Generate a multinomial vector in-place
void DirichletDistribution::generate(Vector& y) const
{
        //Vector y(n);
    real sum = 0.0;
    for (int i=0; i<n; i++) {
        y(i) = gengam(1.0, alpha(i));
        sum += y(i);
    }
    real invsum = 1.0 / sum;
     y *= invsum;
}

/** Dirichlet distribution
    Gets the parameters of a multinomial distribution as input.
*/
real DirichletDistribution::pdf(const Vector& x) const
{
	return exp(log_pdf(x));
}

/** Dirichlet distribution

    Gets the parameters \f$x\f$ of a multinomial distribution as
	input.
	
	
	Returns the logarithm of the pdf.

*/
real DirichletDistribution::log_pdf(const Vector& x) const
{
    assert(x.Size() == n);

    real log_prod = 0.0;
    real sum = 0.0;
    for (int i=0; i<n; i++) {
        real xi = x(i);
        if (xi<=0) {
            Swarning ("Got a negative value for x[%d]:%f\n", i, xi);
            return 0.0;
        }
        sum += xi;
        log_prod += log(xi) * alpha(i);
    }
    if (fabs(sum-1.0f)>0.001) {
        Swarning ("Vector x not a distribution apparently: sum=%f.  Returning 0.\n", sum);
        return 0.0;
    }
    return log_prod - logBeta(alpha);
}

/// Update with a vector of observations sampled from a Multinomial
///
/// Note that this is different from observing a sample from a Multinomial!
void DirichletDistribution::update(Vector* x)
{
	real tmp = 0;
    for (int i=0; i<n; ++i) {
        real xi = (*x)(i);
        alpha(i) += xi;
        alpha_sum += xi;
		tmp += xi;
    }
	n_observations += (int) tmp;
}


/// When there is only one observation, give it directly.
real DirichletDistribution::Observe(int i)
{
    //real p = alpha(i) / alpha.Sum();
    assert(fabs(alpha.Sum() - alpha_sum) < 1e-6);
    real p = alpha(i) / alpha_sum; //alpha.Sum();
    alpha(i) += 1.0;
    alpha_sum += 1.0;
	n_observations++;
    return p;
}

/// Return the parameters
Vector DirichletDistribution::getParameters() const
{
	return alpha;
}

/// Return the marginal probabilities
Vector DirichletDistribution::getMarginal() const
{
    Vector p = alpha;
    if (alpha_sum > 0) {
        p /= alpha_sum;
    } else {
        real invs = 1.0 / (real) p.Size();
        for (int i=0; i<p.Size(); i++) {
            p(i) = invs;
        }
    }
    //printf ("sum: %f\n", p.Sum());
    assert(fabs(p.Sum()-1.0) < 0.0001);
	return p;
}

/// Return the marginal probabilities
real DirichletDistribution::marginal_pdf(int i) const
{
	return alpha(i) / alpha_sum;
}

