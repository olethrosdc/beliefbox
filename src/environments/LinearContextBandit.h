// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTEXT_BANDIT_H
#define CONTEXT_BANDIT_H

#include "DiscreteMDP.h"
#include "Environment.h"
#include "RandomNumberGenerator.h"
#include "NormalDistribution.h"
#include "Matrix.h"
#include <string>
#include <vector>

/** A linear context n-armed bandit.
	
	The reward at time t is \f$r_t \mid a_t = i, x_t = x \sim {\cal N}(x' M_i, \sigma_i)\f$
 */
class LinearContextBandit : public ContinuousStateEnvironment
{
protected:
	std::vector<Vector> mean;
	Vector var;
    Matrix G; ///< generator matrix
    RandomNumberGenerator* rng;
    Vector U_x; ///< upper bound
    Vector L_x; ///< lower bound
public:
    LinearContextBandit(uint n_actions_,
						uint n_features_,
						RandomNumberGenerator* rng_);
    virtual ~LinearContextBandit();
    void GenerateContext();
    const Vector& StateUpperBound() const
    {
        return U_x;
    }
    const Vector& StateLowerBound() const
    {
        return L_x;
    }

    virtual void Reset();
    virtual bool Act(int action);
    virtual const char* Name()
    {
        return "Linear Context Bandit";
    }
protected:
    NormalDistribution normal;
};

#endif
