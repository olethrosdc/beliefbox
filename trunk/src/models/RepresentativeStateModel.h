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

/** RSM 

    Model must implement
    - real getExpectedReward(const S&, const A&);
    - real getTransitionProbability(const S&, const A&);

    In addition, if we do not pass a set of states, we need:
    - const S& StateLowerBound()
    - const S& StateUpperBound()

    Sampler
    - The sampler

    S is the state class
    A is the action class
 */
template <class Model, class S, class A>
class RepresentativeStateModel
{
protected:
    real gamma;					///< discount factor
	real threshold; 			///< Value iteration accuracy threshold
	Model& model;				///< model to use to build discrete approximation from
	uint init_samples;			///< initial number of collected samples 
    std::vector<S> states;		///< set of representative states
	uint n_actions;				///< number of actions (Note: what to do for continuous?)
    Vector V;					///< value vector cache
	DiscreteMDP* mdp;			///< pointer to approximate model
public:
    /// Build a model from a set of representative states
    RepresentativeStateModel(real gamma_, real threshold_, Model& model_, const std::vector<S>& states_, uint n_actions_) :
        gamma(gamma_),
		threshold(threshold_),
        model(model_),
        states(states_),
		n_actions(n_actions_),
        mdp(NULL)
        //mdp(states.size(), n_actions)
    {
		init_samples = states.size();
		assert(n_actions > 0);
        BuildMDP();
		ComputeStateValues();
    }
    /// Build a model from n_states representative states
    RepresentativeStateModel(real gamma_,
							 real threshold_,
							 Model& model_,
                             uint n_states,
                             uint n_actions_) :
        gamma(gamma_),
		threshold(threshold_),
        model(model_),
		n_actions(n_actions_),
        mdp(NULL)
    {
        assert(n_actions > 0);
        S lower_bound = model.StateLowerBound();
        S upper_bound = model.StateUpperBound();
		init_samples = n_states;
        //logmsg("[%d %d]\n", lower_bound, upper_bound);
        for (uint i=0; i<n_states; ++i) {
            //S state = i % (int) upper_bound;
            S state = urandom(lower_bound, upper_bound);
            printf("%d -> %d\n", i, state);
            states.push_back(state);
        }
        BuildMDP();
		ComputeStateValues();
    }


    /// Build a model using a sampler
    template <class Sampler>
    RepresentativeStateModel(real gamma_,
							 real threshold_,
							 Model model_,
                             Sampler& sampler,
                             uint n_states,
                             uint n_actions_) :
        gamma(gamma_),
		threshold(threshold_),
        model(model_),
		n_actions(n_actions_),
        mdp(NULL)
    {
        assert(n_actions > 0);
		init_samples = n_states;
        for (uint i=0; i<n_states; ++i) {
            states.push_back(sampler.Generate());
        }
        BuildMDP();
		ComputeStateValues();
    }

	~RepresentativeStateModel()
	{ 
		if (mdp) {
			delete mdp;
		}
	}
	
	void Reset()
	{
		states.resize(init_samples);
		S lower_bound = model.StateLowerBound();
        S upper_bound = model.StateUpperBound();
		for (uint i=0; i<init_samples; ++i) {
            S state = urandom(lower_bound, upper_bound);
            states.push_back(state);
        }
		BuildMDP();
	}
	
    /// Add a new state to the representative set
    void AddState(const S& state)
    {
        states.push_back(state);
		//BuildMDP();
//		ComputeStateValues();
    }
    
    /// Build the approximate MDP
    void BuildMDP()
    {
        int n_states = states.size();
        mdp = new DiscreteMDP(n_states, n_actions);
        Vector p(n_states);
//        printf("Build model\n");
		for (int i=0; i<n_states; ++i) {
			for (uint a=0; a<n_actions; ++a) {
                // Set rewards
                //logmsg("Getting reward\n");
                real r_ia = model.getExpectedReward(states[i], a);
                //logmsg("Got it!\n");
                mdp->setFixedReward(i, a, r_ia);

                // Set transition probabilities
				for (int j=0; j<n_states; ++j) {
					p(j) =  model.getTransitionProbability(states[i], a, states[j]);
				}
                //logmsg ("s:%d a:%d r:(%f %f) ", i, a, r_ia,  model.getExpectedReward(states[i], a)); p.print(stdout);
                real sum = p.Sum();
				if (sum > 0) {
                    p /= sum;
                } else {
                    //Swarning("OUCH\n");
                    p += 1.0;
                    p /= p.Sum();
                }
				logmsg ("s:%d a:%d r:(%f %f) ", i, a, r_ia,  model.getExpectedReward(states[i], a)); p.print(stdout);
                mdp->setTransitionProbabilities(i, a, p, 1e-3);
			}
		}

        mdp->Check();
    }
	
	void Update(Model& model_)
	{
		model = model_;
		BuildMDP();
		ComputeStateValues();
	}

	void ComputeStateValues(int max_iter = -1)
    {
		ValueIteration value_iteration(mdp, gamma);
		value_iteration.ComputeStateValuesStandard(threshold, max_iter);
        for (int i=0; i<(int) states.size(); ++i) {
            printf("%f\n", value_iteration.getValue(i));
        }
		V = value_iteration.V;
		//printf("Value state\n");
        //logmsg("AV: "); V.print(stdout);
	}

	real getValue(const S& state, const A& action)
	{
		Vector p((int) states.size());
		for (uint j=0; j<states.size(); ++j) {
			p(j) = model.getTransitionProbability(state, action, states[j]);
		}
        real sum = p.Sum();
        
        if (sum > 0) {
            p /= sum;
        } else {
            p += 1.0;
            p /= p.Sum();
        }
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
	
	std::vector<S> getRepresentativesSamples()
	{
		return states;
	}
	
	int getNSamples()
	{
		return states.size();
	}
	S getSample(int i)
	{
		assert(i >= 0 && i < states.size());
		return states[i];
	}
};


#endif
