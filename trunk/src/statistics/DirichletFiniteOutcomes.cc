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

#include "DirichletFiniteOutcomes.h"
#include "ranlib.h"
#include "SpecialFunctions.h"

/// Create a placeholder Dirichlet
DirichletFiniteOutcomes::DirichletFiniteOutcomes()
    : DirichletDistribution(), n_seen_symbols(0)
{
}

/// Create a Dirichlet with uniform parameters
DirichletFiniteOutcomes::DirichletFiniteOutcomes(int n, real p)
    : DirichletDistribution(n, 0.0), prior_alpha(p), n_seen_symbols(0)
{
    alpha_sum = p;
}


/// Destructor
DirichletFiniteOutcomes::~DirichletFiniteOutcomes()
{
}

/// Generate a multinomial vector
Vector DirichletFiniteOutcomes::generate() const
{
    Vector x(n);
    generate(x);
    return x;
}

/// Generate a multinomial vector in-place
void DirichletFiniteOutcomes::generate(Vector& y) const
{
    real N_t = (real) n_seen_symbols;
    //real Z = (1 + N) * prior_alpha + S; // total dirichlet mass

    // generate probability of observing a new symbol
    real ha_0 = gengam(1.0, alpha_sum);
    real ha_1 = gengam(1.0, (1.0 + N_t) * prior_alpha);
    real p_old = ha_0 / (ha_0 + ha_1);
    
    real sum_0 = 0.0;
    real sum_1 = 0.0;
    for (int i=0; i<n; i++) {
        if (alpha(i) > 0) {
            y(i) = gengam(1.0, alpha(i));
            sum_0 += y(i);
        } else {
            y(i) = gengam(1.0, 1.0);
            sum_1 += y(i);
        }
    }
    real is0 = p_old / sum_0;
    real is1 = (1.0 - p_old) / sum_1;

    for (int i=0; i<n; i++) {
        if (alpha(i) > 0) {
            y(i) *= is0;
        } else {
            y(i) *= is1;
        }
    }

    assert(fabs(y.Sum() - 1.0) < 1e-6);

}

/** Dirichlet distribution
    Gets the parameters of a multinomial distribution as input.
*/
real DirichletFiniteOutcomes::pdf(const Vector& x) const
{
    Swarning("Not correctly implemented\n");
	return exp(log_pdf(x));
}

/** Dirichlet distribution

    Gets the parameters \f$x\f$ of a multinomial distribution as
	input.
	
	
	Returns the logarithm of the pdf.

*/
real DirichletFiniteOutcomes::log_pdf(const Vector& x) const
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
    Swarning("Not correctly implemented\n");
    return log_prod - logBeta(alpha);
}

void DirichletFiniteOutcomes::update(Vector* x)
{
    for (int i=0; i<n; ++i) {
        real xi = (*x)(i);
        if (xi > 0 && alpha(i) == 0) {
            n_seen_symbols++;
        }
        alpha(i) += xi;
        alpha_sum += xi;
    }
}

/// When there is only one observation, give it directly.
real DirichletFiniteOutcomes::Observe(int i)
{
    real p = getMarginal()(i);
    if (alpha(i) == 0) {
        n_seen_symbols++;
    }
    alpha(i) += 1.0;
    alpha_sum += 1.0;
    return p;
}


/// Return the marginal probabilities
Vector DirichletFiniteOutcomes::getMarginal() const
{
    real N_t = (real) n_seen_symbols;
    real Z = (1 + N_t) * prior_alpha + alpha_sum; // total dirichlet mass
    Vector P = (alpha + prior_alpha) / Z;

    // distributed the remaining mass equally across unobserved outcomes
    real n_zero_outcomes = n - N_t;
    if (n_zero_outcomes > 0) {
        real SA = 1.0 / n_zero_outcomes;
        for (int i=0; i<n; ++i) {
            if (alpha(i)==0) {
                P(i) *= SA;
            }
        }
    }
    return P;
}

void DirichletFiniteOutcomes::resize(int n, real p)
{
    this->n = n;
    prior_alpha = p;
    alpha.Resize(n);
    alpha.Clear();
}
