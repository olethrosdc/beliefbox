// -*- Mode: c++ -*-
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: PolicyEstimation.c,v 1.5 2006/11/08 17:20:17 cdimitrakakis Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#define _DEBUG_GRADIENT_LEVEL 100

#include "PolicyGradient.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include "GeometricDistribution.h"
#include "MultinomialDistribution.h"
#include <cmath>
#include <cassert>

PolicyGradient::PolicyGradient(const DiscreteMDP* mdp_,
                               real gamma_,
                               real step_size_)
    : evaluation(NULL, mdp_, gamma_),
      mdp(mdp_),
      starting(mdp->getNStates()),
      gamma(gamma_),
      step_size(step_size_)
{
    assert (mdp);
    assert (gamma>=0 && gamma <=1);
    
    n_actions = mdp->getNActions();
    n_states = mdp->getNStates();
    
    policy = new FixedDiscretePolicy(n_states, n_actions);
    real p = 1.0 / (real) n_actions;
    for (int s=0; s<n_states; s++) {
        Vector* Theta = policy->getActionProbabilitiesPtr(s);
        for (int a=0; a<n_actions; a++) {
            (*Theta)(a) = p;
        }
    }

    evaluation.policy = policy;
    
    a_max.resize(n_states);
    for (int i=0; i<n_states; ++i) {
        starting(i) = 1.0 / (real) n_states;
    }
    Reset();
}


void PolicyGradient::Reset()
{
    policy->Reset();
    for (int s=0; s<n_states; s++) {
        a_max[s] = ArgMax(policy->getActionProbabilitiesPtr(s));
    }
}

PolicyGradient::~PolicyGradient()
{
    delete policy;
}

/** Model-based policy gradient
   
    threshold - exit policy gradient is smaller than a threshold
    max_iter - exit when the number of iterations reaches max_iter

*/
void PolicyGradient::ModelBasedGradient(real threshold, int max_iter)
{
    bool policy_stable = true;
    int iter = 0;	
	//Matrix G (n_states, n_actions);
    do {
        policy_stable = false;
        Delta = 0.0;
        // evaluate policy

		evaluation.ComputeStateValues(threshold, max_iter);
		
		Matrix D(n_states, n_actions);
        for (int i=0; i<n_states; i++) {
            real baseline = evaluation.getValue(i);
            for (int a=0; a<n_actions; a++) {
                D(i, a) = starting(i) * (evaluation.getValue(i, a) - baseline);
                Delta += fabs(D(i,a));
            }
        }

        //printf("# D\n"); D.print(stdout);
        //printf(" W\n");
        Delta = 0;
        for (int i=0; i<n_states; i++) {
            Vector* Theta = policy->getActionProbabilitiesPtr(i);
            //Theta->print(stdout);
            real s = 0;
            for (int a=0; a<n_actions; a++) {
                real new_value = (*Theta)(a) + step_size * D(i, a);

                if (new_value < 0) {
                    new_value = 0;
                }
                if (new_value > 1) {
                    new_value = 1;
                }
                s += new_value;
                Delta += fabs((*Theta)(a) - new_value);
                (*Theta)(a) = new_value;
            }
            //printf (" -- %f\n", s);
            (*Theta) /= s;
            //Theta->print(stdout);
        }

        real U = 0;
        for (int i=0; i<n_states; i++) {
            U += starting(i) * evaluation.getValue(i);
        }
        //printf ("%f %f %d # Utility\n", U, Delta, iter);
        if (Delta < threshold && iter > 0) {
            policy_stable = true;
        }

        if (max_iter >= 0) {
            iter++;
        }
    } while(policy_stable == false && iter < max_iter);

    //printf ("delta = %f, iter left = %d\n", Delta, max_iter - iter);
    
    //policy->Show();
}




/** Use the feature expectation version instead.

	The computations here should be equivalent to the once in the other model-based gradient. However, they are inherently more complex, and convergence might be slower.

*/
void PolicyGradient::ModelBasedGradientFeatureExpectation(real threshold, int max_iter)
{
    bool policy_stable = true;
    int iter = 0;	
	//Matrix G (n_states, n_actions);
    do {
        policy_stable = false;
        Delta = 0.0;
        // evaluate policy

        

        evaluation.ComputeStateValuesFeatureExpectation(threshold, max_iter);
                
        Matrix& F = evaluation.FeatureMatrix;
        // r = E[reward | s]
        Vector r(n_states);
        for (int i=0; i<n_states; i++) {
            for (int a=0; a<n_actions; a++) {
                r(i) += mdp->getExpectedReward(i, a) * policy->getActionProbability(i, a);
            }
        }
        // L = y' X
        Vector L(n_states);
        for (int i=0; i<n_states; i++) {
            for (int j=0; j<n_states; j++) {
                L(i) += starting(j) * F(i, j);
            }
        }

        // R = X r
        Vector R(n_states);
        for (int i=0; i<n_states; i++) {
            for (int j=0; j<n_states; j++) {
                R(i) += F(i, j) * r(j);
            }
        }


        Matrix D(n_states, n_actions);
        for (int a=0; a<n_actions; a++) {
            for (int i=0; i<n_states; i++) {
                for (int j=0; j<n_states; j++) {
                    D(i, a) += mdp->getTransitionProbability(i, a, j) * R(j);
                }
                D(i, a) *= L(i);
                Delta += fabs(D(i,a));
            }
        }
		
        //printf(" W\n");
        Delta = 0;
        for (int i=0; i<n_states; i++) {
            Vector* Theta = policy->getActionProbabilitiesPtr(i);
            //Theta->print(stdout);
            real s = 0;
			for (int a=0; a<n_actions; a++) {
                s += D(i, a);
			}
			// The sum of all these should be zero ideally
			real fudge = s / (real) n_actions;;
			s = 0;
            for (int a=0; a<n_actions; a++) {
                real new_value = (*Theta)(a) + step_size * (D(i, a) - fudge);

                if (new_value < 0) {
                    new_value = 0;
                }
                if (new_value > 1) {
                    new_value = 1;
                }
                s += new_value;
                Delta += fabs((*Theta)(a) - new_value);
                (*Theta)(a) = new_value;
            }
            //printf (" -- %f\n", s);
            (*Theta) /= s;
            //Theta->print(stdout);
        }

        real U = 0;
        for (int i=0; i<n_states; i++) {
            U += starting(i) * evaluation.getValue(i);
        }
        //printf ("%f %f %d # Utility\n", U, Delta, iter);
		//policy->Show();

        if (Delta < threshold && iter > 0) {
            policy_stable = true;
        }

        if (max_iter >= 0) {
            iter++;
        }
    } while(policy_stable == false && iter < max_iter);
	
    //printf ("delta = %f, iter left = %d\n", Delta, max_iter - iter);
    //policy->Show();
}


/** Sample trajectories from the model to compute gradients. 
 *
 * This has the advantage that the Markov property is not required, but convergence is worse.
 * Here \f$\nabla U(\pi, \mu) = \sum_h U(h) P_\mu^\pi(h) \sum_t \frac{\nabla \pi(a_t | h_t)}{\pi(a_t | h_t)}\f$.
 */
void PolicyGradient::TrajectoryGradient(real threshold, int max_iter, int n_samples)
{
	MultinomialDistribution starting_state_distribution(starting);
	GeometricDistribution horizon_distribution(1 - gamma);



	bool use_remaining_return = false; 	//< Instead of using U(h), use U(h_{t:T}) for each action.
	bool fast_softmax = true; //< use the fast softmax calculation
	
	// randomly initialise the parameter matrix
	Matrix params(n_states, n_actions); ///< parameters
	for (int s=0; s<n_states; ++s) {
		for (int a=0; a<n_actions; ++a) {
			params(s,a) = urandom(-0.5,0.5);
		}
	}

	Matrix D(n_states, n_actions);
	Matrix D2(n_states, n_actions);
	real Delta = 0;
	for (int iter=0; iter<max_iter; ++iter) {
		D.Clear();
		D2.Clear();
					
		// number of samples before policy evaluation
		real fudge = 0;
		for (int sample=0; sample<n_samples; ++sample) {
			// Get a sample trajectory
			int state = starting_state_distribution.generateInt();
			//int horizon = 1.0f + 1.0f/(1 - gamma); // use this for fixed horizon, but make sure to discount rewards!
			int horizon = horizon_distribution.generate(); // use this for unbiased samples
			std::vector<int> states(horizon);
			std::vector<int> actions(horizon);
			std::vector<real> rewards(horizon);
			real utility = 0;
			policy->Reset(state);
			for (int t=0; t<horizon; t++) {
				int action = policy->SelectAction();
				states[t] = state;
				actions[t] = action;
				real reward = mdp->generateReward(state, action);
				rewards[t] = reward;
				state = mdp->generateState(state, action);
				utility += reward;
			}
			//fudge += utility;
			// calculate the gradient direction

			for (int t=0; t<horizon; t++) {
				if (fast_softmax) {
					// uses the property of the softmax to make a simpler gradient calculation
					for (int a=0; a<n_actions; ++a) {
						real d_sa = utility;
						real p_a = policy->getActionProbability(states[t], a);
						if (a==actions[t]) {
							d_sa *= 1 - p_a;
						} else {
							d_sa *= - p_a;
						}
#if _DEBUG_GRADIENT_LEVEL > 90
						printf("d:%f (s:%d a:%d r:%f u:%f)\n", d_sa, states[t], a, rewards[t], utility);
#endif
						D(states[t], a) += d_sa;
					}
				} else {
					// naive calculation. Should be the same as the
					// one calculating the log directly, but for some reason it's off by a scaling factor for the non-chosen actions.
					Vector eW = exp(params.getRow(states[t]));
					real S = eW.Sum();
					real d_sa = utility / policy->getActionProbability(states[t], actions[t]);
					// softmax - straight-up implementation
					for (int a=0; a<n_actions; ++a) {
						if (a==actions[t]) {
							d_sa *= eW(a)*(S - eW(a)) / (S*S);
						} else {
							d_sa *= - eW(a)*eW(actions[t]) / (S*S);
						}
#if _DEBUG_GRADIENT_LEVEL > 90
						printf("d:%f (s:%d a:%d r:%f u:%f)\n", d_sa, states[t], a, rewards[t], utility);
#endif
						D(states[t], a) += d_sa;
					}
				}
				if (use_remaining_return) {
					utility -= rewards[t];
				}
			}

			// update parameters after multiple trajectory samples
			//D *=1.0 /((real) horizon);

		}
#if _DEBUG_GRADIENT_LEVEL > 0
		printf("---D--- %d/%d---\n", iter, max_iter); D.print(stdout);
		printf("--- params ----\n"); params.print(stdout);
#endif

		// update parameters
		Delta = D.L2Norm() / (real) n_samples;
		params += step_size * (D- fudge) / (real) (n_samples + iter);		
		//printf("eW\n");

#if 0
		// normalise and scale down parameters slightly
		for (int s=0; s<n_states; ++s) {
			real max = Max(params.getRow(s));
			for (int a=0; a<n_actions; ++a) {
				params(s,a) -= max;
			}
		}
		params*=0.999;
#endif
		// update policy from parameters
		for (int s=0; s<n_states; ++s) {
			Vector eW = exp(params.getRow(s));
			eW /= eW.Sum();
			Vector* pS = policy->getActionProbabilitiesPtr(s);
			(*pS) = eW;
		}
		// calculate an evaluation of the policy
		if (1) // evaluate
			{
				evaluation.ComputeStateValuesFeatureExpectation(threshold, max_iter);
				real U = 0;
				for (int i=0; i<n_states; i++) {
					U += starting(i) * evaluation.getValue(i);
				}
				printf ("%f %f %d # Utility\n", U, Delta, iter);
			}
	}
	policy->Show();
}


