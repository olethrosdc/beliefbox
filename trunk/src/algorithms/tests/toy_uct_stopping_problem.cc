// -*- Mode: C++; -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "PolicyEvaluation.h"
#include "BetaDistribution.h"

#include <list>

// We start with an initial belief and then we expand all its possible
// observations.  So, each belief node has to have some
// characteristics.

class BeliefNode
{
public:
    virtual ~BeliefNode()
    {
    }
};

class BetaVectorBeliefNode : public BeliefNode
{
public:
    std::vector<BetaDistribution> posterior = prior;
    virtual ~BetaBeliefNode()
    {
    }
};


/// A toy UCT stopping problem
class BeliefExpansionAlgorithm 
{
public:
    std::vector<BetaDistribution> prior; /// prior for all actions
    uint n_actions;
    BeliefExpansionAlgorithm(std::vector<ConjugatePrior> prior_)
        : prior(prior_)
    {
        n_actions = prior.size();
    }
    virtual ~BeliefExpansionAlgorithm
    {
    }
    void Observe(int action, real reward)
    {
        assert (action>= 0 && action < n_actions);
        real x = (real) (((int) reward) * 2 - 1);
        std::vector<BetaDistribution> posterior = prior;
        posterior[action].calculatePosterior(x);
    }
};

/// A toy UCT stopping problem
class UCTBeliefExpansion : BeliefExpansionAlgorithm
{
public:
    BeliefExpansionAlgorithm()

    virtual ~UCTBeliefExpansion
    {
    }
    int Act()
    {
	
    }
    void Observe(reward);
};



int main (int argc, char** argv)
{
    real alpha = 0;
    real beta = 0;
    ConjugatePrior* prior = new BetaDistribution(alpha, beta);

    real mean_r = urandom(-1, 1);
    
    EvaluateAlgorithm(algorithm, mean_r);
    
    return 0;
}

void EvaluateAlgorithm(UCTAlgorithm* algorithm, real mean_r)
{
    
    // blah
}

#endif
