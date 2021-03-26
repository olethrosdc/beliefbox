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

//#define _LIBCPP_DEBUG_LEVEL 1

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
        printf ("%f %f %d # Utility\n", U, Delta, iter);
		policy->Show();

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
 * This has the advantage that the Markov property is not required.
 * Here \f$\nabla U(\pi, \mu) = \sum_h U(h) P_\mu^\pi(h) \sum_t \frac{\nabla \pi(a_t | h_t)}{\pi(a_t | h_t)}\f$.
 */
void PolicyGradient::TrajectoryGradient(real threshold, int max_iter)
{
	MultinomialDistribution starting_state_distribution(starting);
	GeometricDistribution horizon_distribution(1 - gamma);

	for (int iter=0; iter<max_iter; ++iter) {
		// Get a sample trajectory
		int state = starting_state_distribution.generateInt();
		int horizon = 1 + horizon_distribution.generate(); // use this for unbiased samples
		std::vector<int> states(horizon);
		std::vector<int> actions(horizon);
		real utility = 0;
		policy->Reset(state);
		for (int t=0; t<horizon; t++) {
            int action = policy->SelectAction();
			states[t] = state;
			actions[t] = action;
			real reward = mdp->generateReward(state, action);
			state = mdp->generateState(state, action);
			utility += reward;
		}
		
		// calculate the gradient direction
		Matrix D(n_states, n_actions);
		for (int t=0; t<horizon; t++) {
			Vector* Theta  = policy->getActionProbabilitiesPtr(states[t]);
			real S = (*Theta).Sum();
			for (int a=0; a<n_actions; ++a) {
				real d_sa = utility / (*Theta)(actions[t]);
				if (a==actions[t]) {
					d_sa *= 1; //(S - (*Theta)(a)) / (S*S);
				} else {
					d_sa *= 0; //- (*Theta)(a) / (S*S);
				}
				//printf("%f (%d %d)\n", d_sa, states[t], a);
				Delta += fabs(d_sa);
				D(states[t], a) += d_sa;
			}
		}
		D *=1.0 /((real) horizon);
		
		// Apply the gradient direction
        //printf(" D:\n");
		//D.print(stdout);
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

                if (new_value < 1e-6) {
                    new_value = 1e-6;
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

		
		if (1) // evaluate
		{
			evaluation.ComputeStateValuesFeatureExpectation(threshold, max_iter);
			real U = 0;
			for (int i=0; i<n_states; i++) {
				U += starting(i) * evaluation.getValue(i);
			}
			printf ("%f %f %d # Utility\n", U, Delta, iter);
			//policy->Show();
		}
	}
	policy->Show();

}




