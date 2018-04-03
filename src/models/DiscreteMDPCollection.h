// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_MDP_COLLECTION
#define DISCRETE_MDP_COLLECTION

#include "DiscreteMDPCounts.h"
#include "DiscreteMDPAggregate.h"
#include "Gridworld.h"
#include <set>
#include <vector>

/** A collection of discrete MDPs
    
    
    this is just a collection
 */
class DiscreteMDPCollection : public MDPModel
{
protected:
    std::vector<DiscreteMDPCounts*> A; ///< Pointers to the MDPs
    std::vector<real> P; ///< Probabilities of MDPs
public:
    DiscreteMDPCollection(int n_aggregates, int n_states, int n_actions);
    DiscreteMDPCollection(Gridworld& gridworld, int n_aggregates, int n_states, int n_actions);
	virtual DiscreteMDPCollection* Clone() const
	{
		Serror("Not implemented!\n");
		return NULL;
	}

    virtual ~DiscreteMDPCollection();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual real getExpectedReward (int s, int a) const;
    virtual void Reset();
    int get_n_models();
    std::vector<real> & GetModelProbabilities();
    std::vector<DiscreteMDPCounts*> & GetModels();
    virtual DiscreteMDP* generate() const
    {
        fprintf(stderr, "Not implemented!\n");
        exit(-1);
    }
    virtual const DiscreteMDP* const getMeanMDP() const
    {
        fprintf(stderr, "There is no mean MDP!\n");
        exit(-1);
    }

};



#endif
