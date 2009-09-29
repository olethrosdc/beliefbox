/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteHiddenMarkovModelOnlineEM.h"
#include "RandomNumberGenerator.h"
#include <cassert>
#include <cmath>

DiscreteHiddenMarkovModelOnlineEM::DiscreteHiddenMarkovModelOnlineEM(int n_states_, int n_observations_)  : DiscreteHiddenMarkovModel(n_states, n_observations), q(n_states)
{
    InitPhi();
}


DiscreteHiddenMarkovModelOnlineEM::DiscreteHiddenMarkovModelOnlineEM(Matrix& Pr_S, Matrix& Pr_X) : DiscreteHiddenMarkovModel(Pr_S, Pr_X), q(n_states)
{
    InitPhi();
}

void DiscreteHiddenMarkovModelOnlineEM::InitPhi()
{
    std::vector<int> x(4);
    x[0] = n_states;
    x[1] = n_states;
    x[2] = n_states;
    x[3] = n_observations;
    
    Phi = new Tensor(x);
    q[0] = 1;
}

void DiscreteHiddenMarkovModelOnlineEM::Reset()
{
    DiscreteHiddenMarkovModel::Reset();
    Phi->Reset();
    q.Clear();
    q[0] = 1;
}

DiscreteHiddenMarkovModelOnlineEM::~DiscreteHiddenMarkovModelOnlineEM()
{
        // nothing to do
    delete Phi;
}

/** Update EM estimate
   
    Based on: G. Mongillo and S. Deneve, "Online learning with hidden
    Markov models", Neural computation 20, 1706-1716 (2008).
 */
real DiscreteHiddenMarkovModelOnlineEM::Observe(int x)
{
    assert(x >= 0 && x < n_observations);

    // Calculate gamma
    Matrix gamma(n_states, n_states);
    Matrix AB(n_states, n_states);
    for (int l=0; l<n_states; ++l) {
        for (int h=0; h<n_states; ++h) {
            AB(l,h) = PrS(l, h) * PrX(h, x);
        }
    }
    for (int l=0; l<n_states; ++l) {
        for (int h=0; h<n_states; ++h) {
            real Z = 0;
            for (int m=0; m<n_states; ++m) {
                for (int n=0; n<n_states; ++n) {
                    Z += q[m] * AB(m,n);
                }
            }
            gamma(l, h) = AB(l,h) / Z;
        }
    }

    // update q
    //Matrix gamma_tmp = gamma;
    for (int l=0; l<n_states; ++l) {
        for (int m=0; m<n_states; ++m) {
            //gamma(m, l) = gamma(
        }
    }
    
}

