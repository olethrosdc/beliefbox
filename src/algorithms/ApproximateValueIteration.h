/* -*- Mode: C++; -*- */
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef APPROXIMATE_VALUE_ITERATION_H
#define APPROXIMATE_VALUE_ITERATION_H

#include "ValueFunctionAlgorithm.h"
#include "ValueFunctionModel.h"
#include "Vector.h"

/** Representative state value iteration

    ValueFunctionModel must implement:
    - void AddReturnSample(const S& state, const A& action, const real U); // for inputting a point
	- void CalculateValues();
    - real getValue(const S& state);
    - real getValue(const S& state, const A& action);

    Current Kernel implementations
    -  RBFBasisSet
	
*/
template <typename S, typename A, class VFM, class Model>
class ApproximateValueIteration : public ValueFunctionAlgorithm<S, A>
{
protected:
    real gamma;	///< discount factor
    VFM& vfm; ///< value function model
    std::vector<S> states; ///< representative states
    std::vector<A> actions; ///< representative actions
    Vector V; ///< vector of values
    Model& model; ///< a model
    int n_samples; ///< the number of samples to take to approximate expectations
    int max_iter = 100;
public:
    ApproximateValueIteration(real gamma_,
							  VFM& vfm_,
							  std::vector<S> states_,
							  std::vector<A> actions_,
							  Model& model_,
							  int n_samples_) 
        : gamma(gamma_),
          vfm(vfm_),
          states(states_),
          actions(actions_),
          V((int) states.size()),
          model(model_),
          n_samples(n_samples_)
    {
		max_iter = (int) ceil(1 / (1 - gamma));
    }
    /// Nothing to do in the destrucotr
    virtual ~ApproximateValueIteration()
    {
    }
    virtual void setMaxIter(real max_iter_)
    {
        max_iter = max_iter_;
    }
    /// Specialisation for approximation on a fixed set of states and
    /// actions.
    virtual void CalculateValues()
    {
		// First get an estimate of utilities at the next step
		Matrix Q(states.size(), actions.size());
		for (int iter=0; iter<max_iter; iter++)
            for (uint i=0; i<states.size(); i++) {
                for (uint a=0; a<actions.size(); a++) {
					Q(i, a) = getValue(states.at(i), actions.at(a));
                }
            }

			
		// Fit the samples to the model again
		vfm.Reset();
		for (uint i=0; i<states.size(); i++) {
			for (uint a=0; a<actions.size(); a++) {
				vfm.AddReturnSample(states.at(i),
									actions.at(a),
									Q(i, a));
			}
		}
	}

	
    virtual real getValue(const S& state)  const
    {
		real Q_max = -INF;
		for (uint a=0; a<actions.size(); a++) {
			real Q_a = getValue(state, actions.at(a));
			if (Q_a > Q_max) {
				Q_max = Q_a;
			}
		}
		return Q_max;
    }
	
	// Get the value of a state-action pair usinga finite number of samples
	virtual real getValue(const S& state, const A& action)  const
	{
	real Q_a = 0;
	for (uint k=0; k<n_samples; ++k) {
		model.Reset();
		model.setState(state);
		bool terminal = model.Act(action);
		real r = model.getReward();
		if (terminal) {
			Q_a += r;
		} else {
			S next_state = model.getState();
			// printf ("s: "); state.print(stdout);
			// printf ("r: %f\n", r);
			// printf ("s': "); next_state.print(stdout);
			Q_a += r + gamma * vfm.getValue(next_state);
		}
	}
	return Q_a / (real) n_samples;
}
};

#endif
