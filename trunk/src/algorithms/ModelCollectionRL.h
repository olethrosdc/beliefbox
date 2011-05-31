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

#ifndef MODEL_COLLECTION_RL_H
#define MODEL_COLLECTION_RL_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "ExplorationPolicy.h"
#include "Matrix.h"
#include "real.h"
#include "OnlineAlgorithm.h"
#include "DiscreteMDPCollection.h"
#include "ValueIteration.h"
#include "ModelBasedRL.h"

#include <vector>

class ModelCollectionRL : public ModelBasedRL
{
protected:
    std::vector<const DiscreteMDP*> mdp_vector;
    std::vector<ValueIteration*> vi_vector;
    DiscreteMDPCollection* collection;
    int n_models;
    std::vector<DiscreteMDPCounts*> & models;
public:
    ModelCollectionRL(int n_states_,
                 int n_actions_,
                 real gamma_,
                 real epsilon_,
                 DiscreteMDPCollection* model_,
					  RandomNumberGenerator* rng_,
                 bool use_value_iteration_ = true);
    virtual ~ModelCollectionRL();
    //virtual void Reset();
    /// Full observation
    //virtual real Observe (int state, int action, real reward, int next_state, int next_action);
    /// Partial observation 
    //virtual real Observe (real reward, int next_state, int next_action);
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual int Act(real reward, int next_state);

    virtual real getValue (int state, int action)
    { 
        std::vector<real>& P  = collection->GetModelProbabilities();
        if (use_value_iteration) {
            real Q = 0.0;
            for (int i=0; i<n_models; ++i) {
                Q += P[i] * value_iteration->getValue(state, action);
            }
            return Q;
        } else {
            return 0.0;
        }
    }
    
};

#endif

