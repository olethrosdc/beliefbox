// -*- Mode: c++ -*-
// copyright (c) 2017 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef WEIGHTED_Q_LEARNING_H
#define WEIGHTED_Q_LEARNING_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "Matrix.h"
#include "real.h"
#include "ExplorationPolicy.h"
#include "OnlineAlgorithm.h"
#include "ExponentialDistribution.h"
#include "Distribution.h"
#include <vector>

/** Weighted Ensemble Q-Learning .
	
	This is an online algorithm operating on discrete state-action
	spaces. It is an off-policy algorithm.

        The main parameter is the number of models. For a value equal
        to 1, this is basic Q-learning with a greedy
        policy. Otherwise, online bootstrapping is performed, and a
        version of Thompson sampling is obtained.

	@see Dimitrakakis, 2006: "Ensembles for Sequence Learning"
 */
class WeightedQLearning : public OnlineAlgorithm<int,int>
{
protected:
    int n_models; ///< number of models
    const int n_states; ///< number of states
    const int n_actions; ///< number 
    real gamma; ///< discount factor
    real alpha; ///< learning rate
    real epsilon; ///< probability of selecting a model
    real initial_value; ///< initial value for Q values
    real baseline; ///< baseline reward
    
    std::vector<Matrix> Q; ///< The matrices of Q values

    int current_model; ///< currently used model
    int state; ///< current state
    int action; ///< current action
    int update_interval; ///< update interval
    int t; ///< current time step
    int T; ///< switch horizon
    //ExponentialDistribution weight_distribution; ///< weights distribution for bootstrapping
    BernoulliDistribution weight_distribution; ///< weights distribution for bootstrapping
    Vector weights; ///< bootstrapping weights
public:
    WeightedQLearning(int n_models_,
                      int n_states_,
                      int n_actions_,
                      real gamma_,
                      real alpha_,
                      real epsilon_,
                      real initial_value_= 0.0,
                      real baseline_ = 0.0);
	/// Destructor
    virtual ~WeightedQLearning()
    {
        //ShowModel();
    }
    virtual void Reset();
    virtual real Observe (real reward, int next_state, int next_action);
    virtual real UpdateModel (int m, real reward, int next_state, int next_action);
    virtual int Act(real reward, int next_state);
    void ShowModel()
    {
        for (int i=0; i<n_models; ++i) {
            printf ("model %d\n", i);
            Q[i].print(stdout);
        }
        printf("Current Weights: "); weights.print(stdout);
    }

	/// Get value of state-action
	virtual real getValue (int s, int a)
    {
        real Qsa = 0;
        for (int i=0; i<n_models; ++i) {
            Qsa += Q[i](s, a);
        }
        Qsa /= (real) n_models;
        return Qsa;
    }
};

#endif

