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

#include "DiscreteHiddenMarkovModel.h"

DiscreteHiddenMarkovModel::DiscreteHiddenMarkovModel(int n_states_, int n_observations_)
    : n_states(n_states_), n_observations(n_observations_),
      P_S(n_states), P_X(n_states)
{   
    for (int i=0; i<n_states; ++i) {
        P_S[i].Resize(n_states);
        P_X[i].Resize(n_observations);
    }
    Reset();
}

DiscreteHiddenMarkovModel::DiscreteHiddenMarkovModel(Matrix& Pr_S, Matrix& Pr_X)
    : n_states(Pr_S.Rows()), n_observations(Pr_X.Columns()),
      P_S(n_states), P_X(n_states)
{
    //    printf("# Making HMM with %d states and %d observations\n", n_states, n_observations);
   for (int i=0; i<n_states; ++i) {
        P_S[i].Resize(n_states);
        for (int j=0; j<n_states; ++j) {
            P_S[i].Pr(j) = Pr_S(i, j);
        }
        P_X[i].Resize(n_observations);
        for (int k=0; k<n_observations; ++k) {
            P_X[i].Pr(k) = Pr_X(i, k);
        }
    }
    Reset();
}

void DiscreteHiddenMarkovModel::Reset()
{
    current_state = 0;
}

DiscreteHiddenMarkovModel::~DiscreteHiddenMarkovModel()
{
        // nothing to do
}


int DiscreteHiddenMarkovModel::generate()
{
    current_state = P_S[current_state].generateInt();
    return P_X[current_state].generateInt();
}

int DiscreteHiddenMarkovModel::generate_static()
{
    return P_X[current_state].generateInt();
}

void DiscreteHiddenMarkovModel::Show()
{
  for (int i=0; i<n_states; ++i) {
      for (int j=0; j<n_states; ++j) {
          printf ("%f ", PrS(i, j));
      }
      printf ("# hP_S\n");
  }

  for (int i=0; i<n_states; ++i) {
      for (int j=0; j<n_observations; ++j) {
          printf ("%f ", PrX(i, j));
      }
      printf ("# hP_X\n");
  }
}



DiscreteHiddenMarkovModelStateBelief::DiscreteHiddenMarkovModelStateBelief(DiscreteHiddenMarkovModel* hmm_) : B(hmm_->getNStates()), n_states(hmm_->getNStates()), hmm(hmm_)
{
}

DiscreteHiddenMarkovModelStateBelief::DiscreteHiddenMarkovModelStateBelief(int n_states_) : B(n_states_), n_states(n_states_)
{
}
DiscreteHiddenMarkovModelStateBelief::DiscreteHiddenMarkovModelStateBelief(int n_states_, int start_state) : B(n_states), n_states(n_states_)
{
}
DiscreteHiddenMarkovModelStateBelief::DiscreteHiddenMarkovModelStateBelief(int n_states_, MultinomialDistribution& initial_distribution) : B(initial_distribution), n_states(n_states_)
{
}

/** Update belief from observations
    
    \f[
    b_{t+1}(s_{t+1})
    = b_t(s_{t+1}|x_{t+1})
    = \frac{b_t(x{t+1}|s_{t+1}) b_t(s_{t+1})}{\sum_s b_t(x{t+1}|s_{t+1}=s) b_t(s_{t+1}=s)}
    \f]
*/
real DiscreteHiddenMarkovModelStateBelief::Observe(int x)
{
        //b(s') = sum_i b(s',s=i) = sum_i p(s'|s=i) b(s=i)
    Vector B_next(n_states);
    for (int i=0; i<n_states; ++i) {
        B_next[i] = 0.0;
    }
    for (int src=0; src<n_states; ++src) {
        for (int dst=0; dst<n_states; ++dst) {
            B_next[dst] += hmm->PrS(src, dst) * B.Pr(src);
        }
    }

        //b'(s') = p(x'|s') b(s') / sum_i b(x',s'=i) 
    real sum = 0.0;
    for (int s=0; s<n_states; ++s) {
        B.Pr(s) = hmm->PrX(s, x) * B_next[s];
        sum += B.Pr(s);
    }
    real invsum = 1.0 / sum;
    for (int s=0; s<n_states; ++s) {
        B.Pr(s) *= invsum;
    }
    return sum;
}

real DiscreteHiddenMarkovModelStateBelief::Observe(std::vector<int> x)
{
    real sum = 0.0;
    for (uint i=0; i<x.size(); ++i) {
        sum += Observe(x[i]);
    }
    return sum;
}

Vector DiscreteHiddenMarkovModelStateBelief::getBelief()
{
    return B.getMean();
}

Vector DiscreteHiddenMarkovModelStateBelief::getPrediction()
{
    int n_observations = hmm->getNObservations();
    Vector Px(n_observations);
    for (int x=0; x<n_observations; ++x) {
        Px[x] = 0;
        for (int s=0; s<n_states; ++s) {
            Px[x] += B.Pr(s) * hmm->PrX(s,x);
        }
    }
    return Px;
}

void DiscreteHiddenMarkovModelStateBelief::Reset()
{
    real p = 1.0 / (real) n_states;
    for (int s=0; s<n_states; ++s) {
        B.Pr(s) = p;
    }
}
