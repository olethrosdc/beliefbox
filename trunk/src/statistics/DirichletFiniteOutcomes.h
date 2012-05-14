/* -*- Mode: C++; -*- */
// copyright (c) 20012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DIRICHLET_FINITE_OUTCOMES_H

#include "Dirichlet.h"

/** Finite outcome Dirichlet.

    This is a Dirichlet distirbution with a finite, but unknown, set
    of outcomes. It is a straightforward extension of the Dirichlet
    distribution.

 */
class DirichletFiniteOutcomes : public DirichletDistribution
{
protected:
    real prior_alpha;
    real alpha_sum;
    int n_seen_symbols;
  public:
    DirichletFiniteOutcomes();
    DirichletFiniteOutcomes(int n, real p = 1.0);
    virtual ~DirichletFiniteOutcomes();
    virtual void generate(Vector& x) const;
    virtual Vector generate() const;
    virtual real pdf(const Vector& x) const;
    virtual real log_pdf(const Vector& x) const;
    virtual void update(Vector* x);
    virtual real Observe(int i);
    virtual Vector getMarginal() const;
    virtual void resize(int n, real p = 0.0);

};

#endif


