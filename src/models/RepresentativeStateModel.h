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
    real gamma;
    const Model& model;
    std::vector<S> states;
	uint n_actions;
    Vector V;
	DiscreteMDP* mdp;
public:
    /// Build a model from a set of representative states
    RepresentativeStateModel(real gamma_, const Model& model_, const std::vector<S>& states_, uint n_actions_) :
        gamma(gamma_),
        model(model_),
        states(states_),
		n_actions(n_actions_),
        mdp(states.size(), n_actions)
    {
		assert(n_actions > 0);
        BuildMDP();
    }
    /// Build a model from n_states representative states
    RepresentativeStateModel(real gamma_,
                             const Model& model_,
                             uint n_states,
                             uint n_actions_) :
        gamma(gamma_),
        model(model_),
		n_actions(n_actions_)
    {
        assert(n_actions > 0);
        S lower_bound = model.StateLowerBound();
        S upper_bound = model.StateUpperBound();
        logmsg("[%d %d]\n", lower_bound, upper_bound);
        for (uint i=0; i<n_states; ++i) {
            S state = urandom(lower_bound, upper_bound);
            logmsg("%d -> %d\n", state, i);
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
                // Set rewards
                real r_ia = model.getExpectedReward(states[i], a);
                mdp->setFixedReward(i, a, r_ia);

                // Set transition probabilities
				Vector p(n_states);
				for (int j=0; j<n_states; ++j) {
					p(j) = model.getTransitionProbability(states[i], a, states[j]);
				}
                logmsg ("s:%d a:%d r:(%f %f) ", i, a, r_ia,  model.getExpectedReward(states[i], a)); p.print(stdout);
                real sum = p.Sum();
				if (sum > 0) {
                    p /= sum;
                } else {
                    Swarning("OUCH\n");
                    p += 1.0;
                    p /= p.Sum();
                }

				for (int j=0; j<n_states; ++j) {
					mdp->setTransitionProbability(i, a, j, p(j));
				}
			}
		}
        mdp->Check();
    }

	void ComputeStateValues(real threshold, int max_iter = -1)
    {
		ValueIteration value_iteration(mdp, gamma);
		value_iteration.ComputeStateValues(threshold, max_iter);
		V = value_iteration.V;
        logmsg("AV: "); V.print(stdout);
	}

	real getValue(const S& state, const A& action)
	{
		Vector p((int) states.size());
		for (uint j=0; j<states.size(); ++j) {
			p(j) = model.getTransitionProbability(state, action, states[j]);
		}
		p /= p.Sum();
        real r = model.getExpectedReward(state, action);
        real U = Product(p, V);
		return  r + gamma * U;
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
