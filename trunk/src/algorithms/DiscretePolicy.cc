/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscretePolicy.h"
//#include "MDPDistribution.h"
#include "DiscreteMDP.h"

#include "Random.h"
#include <cmath>

/// Create an uniform discrete policy with a given number of states
/// and actions.
FixedDiscretePolicy::FixedDiscretePolicy(int n_states, int n_actions)
    : DiscretePolicy(n_states, n_actions)
{
    state = 0;
    real p_u = 1.0 / (real) n_actions;
    p.resize(n_states);
    for (uint i=0; i<p.size(); i++) {
        Vector& V = p[i];
        V.Resize(n_actions);
        assert (V.Size()==n_actions);
        for (int j=0; j<n_actions; j++) {
            V[j] = p_u;
        }
    }
}

/// Create a fixed discrete policy from a demonstrator policy
FixedDiscretePolicy::FixedDiscretePolicy(int n_states, int n_actions, Demonstrations<int, int>& D)
    : DiscretePolicy(n_states, n_actions)      
{

    p.resize(n_states);
    for (uint i=0; i<p.size(); i++) {
        p[i].Resize(n_actions); /////!!!!?????
    }

    for (int s=0; s<n_states; ++s) {
        p[s].Resize(n_actions);
        for (int a=0; a<n_actions; ++a) {
            p[s](a) = 1.0 / n_actions;
        }
    }
    for (uint k=0; k<D.trajectories.size(); ++k) {
        for (uint t=0; t<D.trajectories[k].x.size(); ++t) {
            int s = D.trajectories[k].x[t].first;
            int a = D.trajectories[k].x[t].second;
            assert(s >= 0 && s < n_states);
            assert(a >= 0 && a < n_actions);
            p[s](a) += 1.0;
        }
    }

    for (int s=0; s<n_states; ++s) {
        p[s] /= p[s].Sum();
    }

}


/// Create a fixed discrete policy from a set of probability vectors.
FixedDiscretePolicy::FixedDiscretePolicy (std::vector<Vector>& p)
    : DiscretePolicy(n_states, n_actions)
{
    state = 0;
    this->p = p;
    for (uint i=0; i<p.size(); i++) {
        assert (fabs(this->p[i].Sum() - 1.0) < 0.00001);
    }
}


/// Create a greedy policy from a Q value matrix.
FixedDiscretePolicy::FixedDiscretePolicy (int n_states, int n_actions,
                                          Matrix& Q)
    : DiscretePolicy(n_states, n_actions)
{
    assert(Q.Rows() == n_states);
    assert(Q.Columns() == n_actions);

    p.resize(n_states);
    for (uint i=0; i<p.size(); i++) {
        p[i].Resize(n_actions);
    }

    
    for (int s=0; s<n_states; s++) {
        real max_Qa = Q(s, 0);
        int argmax_Qa = 0;
        for (int a=1; a<n_actions; a++) {
            real Qa = Q(s, a);
            if (Qa > max_Qa) {
                max_Qa = Qa;
                argmax_Qa = a;
            }
        }
        Vector* p = getActionProbabilitiesPtr(s);
        for (int a=0; a<n_actions; a++) { 
            (*p)(a) = 0.0;
        }
        (*p)(argmax_Qa) = 1.0;
    }
 }

/// Destructor. Does nothing at the moment.
FixedDiscretePolicy::~FixedDiscretePolicy()
{
}


FixedDiscretePolicy FixedDiscretePolicy::MakeGreedyPolicy()
{
    int n_states = p.size();
    int n_actions = p[0].Size();
    Matrix Q(n_states, n_actions);
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Q(s,a) =getActionProbability(s,a);
        }
    }
    return FixedDiscretePolicy(n_states, n_actions, Q);
}
int FixedDiscretePolicy::SelectAction()
{
    
    int n = p[state].Size();
    real x = urandom();
    real s = 0.0;
    for (int a=0; a<n; ++a) {
        s += p[state][a];
        if (s>x) {
            return a;
        }
    }
    return n-1;
}

void FixedDiscretePolicy::Observe (int& previous_state, int& action, real r, int& next_state)
{
    state = next_state;
}

void FixedDiscretePolicy::Observe (real r, int& next_state)
{
    state = next_state;
}


void FixedDiscretePolicy::Reset(int& start_state)
{
    state = start_state;
}

real FixedDiscretePolicy::getActionProbability(int& action) const
{
    assert(state >= 0 && state < n_states);
    assert(action >= 0 && action < n_actions);
	return p[state][action];
}


real FixedDiscretePolicy::getActionProbability(int& state, int& action) const
{
    assert(state >= 0 && state < n_states);
    assert(action >= 0 && action < n_actions);
	return p[state][action];
}

void FixedDiscretePolicy::Show()
{
    for (int i=0; i != (int) p.size(); ++i) {
        //printf ("%d ", i);
        for (int j=0; j<p[i].Size(); ++j) {
            printf ("%f ", getActionProbability(i, j));
        }
        printf("\n");
    }
}


FixedSoftmaxPolicy::FixedSoftmaxPolicy (Matrix& Q, real beta)
    : FixedDiscretePolicy(Q.Rows(), Q.Columns())
{
    
    for (int s=0; s<Q.Rows(); s++) {
        Vector* p = getActionProbabilitiesPtr(s);
        Vector Qs = Q.getRow(s);
        SoftMax(Qs, *p, beta);
    }
 }

FixedSoftmaxPolicy::~FixedSoftmaxPolicy()
{
}

/** Return a discounted state occupancy Matrix

    Return matrix Phi with
    \f[
    \phi_{i,j} = \sum_{t=0}^T \gamma^t P(s_{t} = j \mid s_0 = i).
    \f]
 */
Matrix DiscountedStateOccupancy(const DiscreteMDP& mdp,
                                const FixedDiscretePolicy& policy,
                                real gamma, real epsilon)
{
    int n_states = mdp.getNStates();
    int n_actions = mdp.getNActions();
    // First calculate the state-state probabilities
    Matrix P(n_states, n_states); 
    P.Clear();
    for (int i=0; i<n_states; ++i) {
        Vector Pi = policy.getActionProbabilities(i);
        for (int a=0; a<n_actions; ++a) {
            real P_a = Pi(a);
            for (int j=0; j<n_states; ++j) {
                real p = mdp.getTransitionProbability(i, a, j);
                P(j, i) += P_a * p;
            }
        }
    }

    logmsg ("Markov chain:\n");
    P.print(stdout);

    const Matrix& Pc = P;    
    // Calculate mu If D is a Nx1 matrix, P is a N*N matrix, then PD
    // is N x 1. We start with D being a singular distribution.
    int T = 1000 + (int) (2.0 * log(epsilon * (1.0 - gamma)) / log(gamma));
    logmsg("Setting horizon to: %d\n", T);
    Matrix Occupancy(n_states, n_states);
    for (int s=0; s<n_states; ++s) {
        Vector mu(n_states);
        real discount = 1;
        Vector D(n_states);
        D(s) = 1.0;
        for (int t=0; t<T; ++t) {
            mu += D * discount;
            const Vector& Dc = D;
            D = Pc * Dc;
            discount *= gamma;
        }
        Occupancy.setRow(s, mu);
    }
    
    return Occupancy;
}

