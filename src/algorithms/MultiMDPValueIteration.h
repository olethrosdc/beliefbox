// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MULTI_MDP_VALUE_ITERATION_H
#define MULTI_MDP_VALUE_ITERATION_H

#include "DiscreteMDP.h"
#include "DiscretePolicy.h"
#include "real.h"
#include <vector>
#include "Matrix.h"

class MultiMDPValueIteration
{
public:
    Vector w;
    std::vector<const DiscreteMDP*> mdp_list;
    real gamma;
    int n_states;
    int n_actions;
    Vector V; ///< actual value
    Vector dV;///< actual delta value
    Vector pV; ///< previous value
    Matrix Q; ///< actual q-value
    Matrix dQ; ///< delta q-value
    Matrix pQ; ///< previous q-value
    int n_mdps; ///< The number of MDPs
    std::vector<Vector> V_mu; ///< the MDPs individual value functions
    real Delta;
    MultiMDPValueIteration(const Vector& w,
                           const std::vector<const DiscreteMDP*>& mdp_list_,
                           real gamma_);
    ~MultiMDPValueIteration();
    void Reset();

    void ComputeStateValues(real threshold, int max_iter=-1);
    void ComputeStateActionValues(real threshold, int max_iter=-1);
    inline real getValue (int state, int action)
    {
        assert(state>=0 && state < n_states);
        assert(action>=0 && action < n_actions);
        return Q(state,action);
    }
    inline real getValue (int state)
    {
        assert(state>=0 && state < n_states);
        return V(state);
    }
    FixedDiscretePolicy* getPolicy();
    void setMDPList(const std::vector<const DiscreteMDP*>& mdp_list_)
    {
        mdp_list = mdp_list_;

        assert(mdp_list.size() == n_mdps);
        assert(w.Size() == n_mdps);

        real w_i = 1.0 / (real) n_mdps;
        for (int i=0; i<n_mdps; i++) {
            w(i) = w_i;
        }
    }
    void SetMDPList(const Vector& w_,
                    const std::vector<const DiscreteMDP*>& mdp_list_)
    {
        w = w_;
        mdp_list = mdp_list_;

        assert(mdp_list.size() == n_mdps);
        assert(w.Size() == n_mdps);
    }
protected:
    real ComputeActionValueForMDPs(int s, int a);
    real ComputeStateActionValueForSingleMDP(int mu, int s, int a);

};

#endif


