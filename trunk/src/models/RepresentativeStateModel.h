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

template <class Model, class S, class A>
class RepresentativeStateModel : public MDP<S, A>
{
protected:
    Model model;
    std::vector<S> states;
	uint n_actions;
    Vector V;
	DiscreteMDP* mdp;

public:
    RepresentativeStateModel(const Model& model_, const std::vector<S>& states_, uint n_actions_) :
        model(model_),
        states(states_),
		n_actions(n_actions_),
		mdp(NULL)
    {
		assert(n_actions > 0);
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
    
	void ComputeStateValues(real gamma, real threshold, int max_iter = -1)
    {
		if (mdp) {
			delete mdp;
		}

        mdp = new DiscreteMDP (states.size(), n_actions);
	
		for (uint i=0; i<states.size(); ++i) {
			for (uint a=0; a<n_actions; ++a) {
				Vector p(states.size());
				for (uint j=0; j<states.size(); ++j) {
					p(j) = model.getTransitionProbability(i, a, j);
				}
				p /= p.Sum();
				for (uint j=0; j<states.size(); ++j) {
					mdp->setTransitionProbability(i, a, j, p(j));
				}
			}
		}

		ValueIteration value_iteration(mdp, gamma);
		value_iteration.ComputeStateValues(threshold, max_iter);
		V = value_iteration.V;
	}

	real getValue(const S& state, const A& action)
	{
		Vector p(states.size());
		for (uint j=0; j<states.size(); ++j) {
			p(j) = model.getTransitionProbability(state, action, j);
		}
		p /= p.Sum();
		return p * V;
	}
    
};


#endif
