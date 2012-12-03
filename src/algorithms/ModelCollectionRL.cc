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

#include "ModelCollectionRL.h"

ModelCollectionRL::ModelCollectionRL(int n_states_,
                                     int n_actions_,
                                     real gamma_,
                                     real epsilon_,
                                     DiscreteMDPCollection* model_,
									 RandomNumberGenerator* rng_,
                                     bool use_value_iteration_) :
    ModelBasedRL(n_states_, n_actions_, gamma_, epsilon_, model_, rng, use_value_iteration_),
    collection(model_),
    n_models(collection->get_n_models()),
    models(collection->GetModels())
{
    mdp_vector.resize(n_models);
    vi_vector.resize(n_models);
    for (int i = 0; i < n_models; ++i) {
        mdp_vector[i] = models[i]->getMeanMDP();
        vi_vector[i] = new ValueIteration(mdp_vector[i], gamma);
    }
}

ModelCollectionRL::~ModelCollectionRL()
{
    for (int i = 0; i < n_models; ++i) {
        delete vi_vector[i];
        //delete mdp_vector[i]; // no deletion! this is a constant!
    }
}

//void ModelCollectionRL::Reset()
//{
//    state = -1;
//    model->Reset();
//}

/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int ModelCollectionRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);

    // update the model
    if (state >= 0) {
        model->AddTransition(state, action, reward, next_state);
    }
    state = next_state;

    std::vector<real>& P  = collection->GetModelProbabilities();
    std::vector<DiscreteMDPCounts*>& M  = collection->GetModels();
    for (int a = 0; a < n_actions; ++a) {
        tmpQ[a] = 0.0;
    }
    
    for (int i = 0; i < n_models; ++i) {
        
        mdp_vector[i] = M[i]->getMeanMDP();

        vi_vector[i]->setMDP(mdp_vector[i]);
        vi_vector[i]->ComputeStateActionValues(0.00, 1);
        for (int a=0; a<n_actions; a++) {
            tmpQ[a] += P[i] * value_iteration->getValue(next_state, a);
        }
    }

    // choose action
    int next_action;
    if (urandom()<epsilon) {
        next_action = (int) floor(urandom(0.0, (real) n_actions));
        //printf ("\n");
    } else {
        next_action = ArgMax(tmpQ);
        //printf (" * A = %d\n", next_action);
    }
    action = next_action;

    return action;
}


