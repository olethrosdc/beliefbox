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
#include "RandomNumberGenerator.h"
#include <cassert>
#include <cmath>

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

real DiscreteHiddenMarkovModel::Expectation(std::vector<int>& observations, Matrix& forward_belief, Matrix& backward_belief)
{
    int T = observations.size();
    assert (forward_belief.Rows() == T);
    assert (forward_belief.Columns() == n_states);

    for (int t=0; t<T; ++t) {
        for (int i=0; i<n_states; ++i) {
            forward_belief(t, i) = -1;
            backward_belief(t, i) = -1;
        }
    }
    //Show();
    
    // initialise the first belief
    Vector B_prev(n_states);
    B_prev[0]= 1.0;
        
    // calculate forward pass
    for (int t=0; t<T; ++t) {
        //b(s') = sum_i b(s',s=i) = sum_i p(s'|s=i) b(s=i)
        Vector B_next(n_states);
        if (t==0) {
            for (int src=0; src<n_states; ++src) {
                for (int dst=0; dst<n_states; ++dst) {
                    B_next[dst] += PrS(src, dst) * B_prev[src];
                }
            }
        } else {    
            for (int src=0; src<n_states; ++src) {
                for (int dst=0; dst<n_states; ++dst) {
                    B_next[dst] += PrS(src, dst) * forward_belief(t - 1, src);
                }
            }
        }

        //b'(s') = p(x'|s') b(s') / sum_i b(x',s'=i) 
        real sum = 0.0;
        for (int s=0; s<n_states; ++s) {
            forward_belief(t, s) = PrX(s, observations[t]) * B_next[s];
            sum += forward_belief(t,s);
        }
        assert (!isnan(sum));
        assert (sum > 0);
        real invsum = 1.0 / sum;
        for (int s=0; s<n_states; ++s) {
            forward_belief(t, s) *= invsum;
        }
        //B_next.print(stdout);
    }
    
    //printf ("Expectation forward pass\n");
    //forward_belief.print(stdout);

    for (int s=0; s<n_states; ++s) {
        backward_belief(T - 1, s) = forward_belief(T - 1, s);
    }

    // calculate backward pass
    Matrix MB(n_states, n_states);
    for (int t=T-2; t>=0; --t) {

        // P(s|x^T) = sum_j P(s,s'=j | x^T)
        //          = sum_j P(s|s'=j,x^T) P(s'=j|x^T)
        for (int j=0; j<n_states; ++j) {
            real sum = 0.0;
            for (int k=0; k<n_states; ++k) {
                MB(k, j) = PrS(k, j) * forward_belief(t, k);
                sum += MB(k, j);
            }
            if (sum == 0) {
                real inv = 1.0 / (real) n_states;
                for (int k=0; k<n_states; ++k) {
                    MB(k, j) = inv;
                }
            } else {                
                real invsum = 1.0 / sum;
                for (int k=0; k<n_states; ++k) {
                    MB(k, j) *= invsum;
                }
            }
        }
        
        // P(s|x^T) = sum_j P(s|x_T) P(s'=j|x^T)
        real sum = 0.0;
        for (int k=0; k<n_states; ++k) {
            backward_belief(t, k) = 0;
            for (int j=0; j<n_states; ++j) {
                backward_belief(t, k) += MB(k, j) * backward_belief(t + 1, j);
            }
            sum += backward_belief(t, k);
        }
    }

    real log_likelihood = 0;
    for (int t=0; t<T; ++t) {
        real p = 0;
        for (int s=0; s<n_states; ++s) {
            p += backward_belief(t, s) * PrX(s, observations[t]);
        }
        log_likelihood += log(p);
    }

    //printf("log likelihood: %f\n", log_likelihood);
    //printf ("Expectation backward pass\n");
    //belief.print(stdout);
    return log_likelihood;
}

void DiscreteHiddenMarkovModel::Maximisation(Matrix& forward_belief, Matrix& backward_belief)
{
    // states first
    // a_ij = sum_t P(x_{t-1} = i, x_t = j | y)/sum_t P(x_{t-1} = i | y)
    int T = forward_belief.Rows();
    assert (forward_belief.Columns() == n_states);

    Matrix hPS(n_states, n_states);

    for (int i=0; i<n_states; ++i) {
        for (int j=0; j<n_states; ++j) {
            real P_ss = 0;
            real P_s = 0;
            for (int t=0; t < T - 1; ++t) {
                real Z = 0;
                for (int k=0; k<n_states; ++k ) {
                    Z += forward_belief(t, k) * PrS(k, j);
                }
                P_ss += forward_belief(t, i) * PrS(i, j) * backward_belief(t + 1, j) / Z;
                P_s += backward_belief(t, i);
            }
            hPS(i,j) = P_ss / P_s;
        }
    }
    
    // b_jk = sum_t P(x_t = j, y_t = k | y) / sum_t P(x_t = j | y) 
    Matrix hPX(n_states, n_observations);
    for (int j=0; j<n_states; ++j) {
        for (int k=0; k<n_observations; ++k) {
            real P_sx = 0;
            real P_s = 0;
            for (int t=0; t < T; t++) {
                P_sx += backward_belief(t, j) * PrX(j, k);
                P_s += backward_belief(t,j);
            }
            hPX(j, k) = P_sx / P_s;
        }
    }

    //hPS.print(stdout);
    //hPX.print(stdout);

    // copy results
    for (int i=0; i<n_states; ++i) {
        assert(approx_eq(hPS.RowSum(i), 1));
        assert(approx_eq(hPX.RowSum(i), 1));
        for (int j=0; j<n_states; ++j) {
            PrS(i, j) = hPS(i, j);
        }
        for (int j=0; j<n_observations; ++j) {
            PrX(i, j) = hPX(i, j);
        }
    }
}


/** Expectation maximisation for HMMs
 
    Calculate
    \f[
    P(x^t|y^t) 
    \f]
    
 */
void DiscreteHiddenMarkovModel::ExpectationMaximisation(std::vector<int>& observations, int n_iterations)
{

    int T = observations.size();
    _belief.Resize(T, n_states);
    Matrix forward_belief(T, n_states);

    for (int iter=0; iter<n_iterations; ++iter) {
        // Expectation step.
        Expectation(observations, forward_belief, _belief);

        // maximisation step
        Maximisation(forward_belief, _belief);
    }
}


//----------------------------------------------------------------------//

DiscreteHiddenMarkovModel* MakeRandomDiscreteHMM(int n_states, int n_observations, real stationarity, RandomNumberGenerator* rng)
{
    assert (n_states > 0);
    assert (n_observations > 0);
    assert (stationarity >= 0 && stationarity <= 1);

    Matrix Pr_S(n_states, n_states);
    Matrix Pr_X(n_states, n_observations);
    for (int i=0; i<n_states; ++i) {
        real sum = 0.0;
        for (int j=0; j<n_observations; ++j) {
            Pr_X(i,j) = 0.1*rng->uniform();
            if (i==j) {
                Pr_X(i,j) += 1.0;
            }
            sum += Pr_X(i,j);
        }
        for (int j=0; j<n_observations; ++j) {
            Pr_X(i,j) /=  sum;
        }
        
    }
    Matrix S(n_states, n_states);
    for (int src=0; src<n_states; ++src) {
        Vector P(n_states);
        for (int i=0; i<n_states; ++i) {
            if (i<=src) {
                P[i] = exp(rng->uniform());
            } else {
                P[i] = exp(10.0*rng->uniform());
            }
        }
        P /= P.Sum();
        P *= (1 - stationarity);
        P[src] += stationarity;
        P /= P.Sum();
            //real sum = 0.0;
        for (int dst=0; dst<n_states; ++dst) {
            Pr_S(src,dst) = P[dst];
        }
    }

    return new DiscreteHiddenMarkovModel (Pr_S, Pr_X);
}

//----------------------------------------------------------------------//

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

/** Update belief from next observation
    
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

/** Update belief from a sequence of observations
 */
real DiscreteHiddenMarkovModelStateBelief::Observe(std::vector<int>& x)
{
    real sum = 0.0;
    for (uint i=0; i<x.size(); ++i) {
        sum += Observe(x[i]);
    }
    return sum;
}

/// Get the current belief
Vector DiscreteHiddenMarkovModelStateBelief::getBelief()
{
    return B.getMean();
}

/// Get the current prediction
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

/// Reset the belief to uniform
void DiscreteHiddenMarkovModelStateBelief::Reset()
{
    real p = 1.0 / (real) n_states;
    for (int s=0; s<n_states; ++s) {
        B.Pr(s) = p;
    }
}
