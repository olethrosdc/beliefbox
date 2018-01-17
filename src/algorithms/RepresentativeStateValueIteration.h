/* -*- Mode: C++; -*- */
// copyright (c) 2012 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
// copyright (c) 2018 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REPRESENTATIVE_STATE_VALUE_ITERATION_H
#define REPRESENTATIVE_STATE_VALUE_ITERATION_H

#include "ValueFunctionAlgorithm.h"
#include "Vector.h"

/** Representative state value iteration

    Kernel must implement:
    - void Evaluate(const S& x); // for inputting a point
    - S F(); // for getting the features corresponding to a point

    Current Kernel implementations
    -  RBFBasisSet
	
*/
template <typename S, typename A, class Kernel, class Model>
class RepresentativeStateValueIteration : public ValueFunctionAlgorithm<S, A>
{
protected:
    real gamma;	///< discount factor
    Kernel& kernel; ///< similarity kernel
    std::vector<S> states; ///< representative states
    std::vector<A> actions; ///< representative actions
    Vector V; ///< vector of values
    Model& model; ///< a model
    int n_samples; ///< the number of samples to take to approximate expectations
    real threshold = 1e-6;
    int max_iter = -1;
public:
    RepresentativeStateValueIteration(real gamma_,
                                      Kernel& kernel_,
                                      std::vector<S> states_,
                                      std::vector<A> actions_,
                                      Model& model_,
                                      int n_samples_) 
        : gamma(gamma_),
          kernel(kernel_),
          states(states_),
          actions(actions_),
          V((int) states.size()),
          model(model_),
          n_samples(n_samples_)
    {
    }
    /// Nothing to do in the destrucotr
    virtual ~RepresentativeStateValueIteration()
    {
    }
    virtual void setThreshold(real threshold_)
    {
        threshold = threshold_;
    }
    virtual void setMaxIter(int max_iter_)
    {
        max_iter = max_iter_;
    }
    /// Specialisation for discrete actions, arbitrary states
    virtual void CalculateValues()
    {
        Vector Vn((int) states.size());
        model.Reset();
        int iter = 0;
        while (1) {
            for (uint i=0; i<states.size(); i++) {
                real Q_max = -INF;
                for (uint a=0; a<actions.size(); a++) {
                    real Q_a = getValue(states.at(i), actions.at(a));
                    //printf ("Q[%d] = %f\n", a, Q_a);
                    if (Q_a > Q_max) {
                        Q_max = Q_a;
                    }
                }
                Vn(i) = Q_max;
            }
            real error = (V - Vn).L1Norm();
            V = Vn;
            //printf ("%d %f\n", iter, error);
            //V.print(stdout);
            if ((max_iter >= 0 && ++iter >= max_iter) || error < threshold) {
                logmsg ("Exiting after %d iterations with %f delta\n", iter, error);
                break;
            }
        }
    }

    real getValue(const S& state)  const
    {
        kernel.Evaluate(state);
        Vector F = kernel.F();
        return Product(V , F / F.Sum());
    }

    /// Get the value of a state-action pair usinga finite number of samples
    real getValue(const S& state, const A& action)  const
    {
        real Q_a = 0;
        for (uint k=0; k<n_samples; ++k) {
            model.Reset();
            model.setState(state);
            bool action_ok = model.Act(action);
            real r = model.getReward();
			if (!action_ok) {
				Q_a += r;
			} else {
				S next_state = model.getState();
				// printf ("s: "); state.print(stdout);
				// printf ("r: %f\n", r);
				// printf ("s': "); next_state.print(stdout);
				Q_a += r + gamma * getValue(next_state);
			}
        }
        return Q_a / (real) n_samples;
    }
};

#endif
