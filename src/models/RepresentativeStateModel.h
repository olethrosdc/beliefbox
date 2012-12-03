/// -*- Mode: c++ -*-
// copyright (c) 2012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REPRESENTATIVE_STATE_MODEL_H
#define REPRESENTATIVE_STATE_MODEL_H

#include "MDP.h"
#include "DiscreteMDP.h"
#include "ValueIteration.h"
#include "Vector.h"


template <class Model, class S, class A>
class RepresentativeStateModel
{
protected:
    const Model& model;
    std::vector<S> states;
	uint n_actions;
    Vector V;
	DiscreteMDP* mdp;
public:
    /// Build a model from a set of representative states
    RepresentativeStateModel(const Model& model_, const std::vector<S>& states_, uint n_actions_) :
        model(model_),
        states(states_),
		n_actions(n_actions_),
        mdp(states.size(), n_actions)
    {
		assert(n_actions > 0);
        BuildMDP();
    }
    /// Build a model from n_states representative states
    RepresentativeStateModel(const Model& model_,
                             uint n_states,
                             uint n_actions_) :
        model(model_),
		n_actions(n_actions_)
    {
        assert(n_actions > 0);
        S lower_bound = model.StateLowerBound();
        S upper_bound = model.StateUpperBound();
        for (uint i=0; i<n_states; ++i) {
            S state = urandom(lower_bound, upper_bound);
            states.push_back(state);
        }
        BuildMDP();
    }

	~RepresentativeStateModel()
	{ 
		if (mdp) {
			delete mdp;
		}
	}
	
    void AddState(const S& state)
    {
        states.push_back(state);
    }
    
    void BuildMDP()
    {
        int n_states = states.size();
        mdp = new DiscreteMDP (n_states, n_actions);
	
		for (int i=0; i<n_states; ++i) {
			for (uint a=0; a<n_actions; ++a) {
				Vector p(n_states);
				for (int j=0; j<n_states; ++j) {
					p(j) = model.getTransitionProbability(i, a, j);
				}
				p /= p.Sum();
				for (int j=0; j<n_states; ++j) {
					mdp->setTransitionProbability(i, a, j, p(j));
				}
			}
		}
    }

	void ComputeStateValues(real gamma, real threshold, int max_iter = -1)
    {
		ValueIteration value_iteration(mdp, gamma);
		value_iteration.ComputeStateValues(threshold, max_iter);
		V = value_iteration.V;
	}

	real getValue(const S& state, const A& action)
	{
		Vector p((int) states.size());
		for (uint j=0; j<states.size(); ++j) {
			p(j) = model.getTransitionProbability(state, action, j);
		}
		p /= p.Sum();
		return Product(p, V);
	}

    
    real getValue(const S& state)
	{
        Vector Q(n_actions);
        for (uint i=0; i<n_actions; ++i) {
            Q(i) = getValue(state, i);
        }
        return Max(Q);
	}

};


#endif
