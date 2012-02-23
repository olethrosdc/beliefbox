/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.c,v 1.6 2006/10/23 16:48:41 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteBanditPolicy.h"
#include "Random.h"

/// Make e-greedy
EpsilonGreedyPolicy::EpsilonGreedyPolicy(int n_actions, real epsilon, ActionValueEstimate* estimator)
{
    this->n_actions = n_actions;
    this->epsilon = epsilon;
    this->estimator = estimator;
}

/// Reset
void EpsilonGreedyPolicy::Reset() {
    estimator->Reset();
}

/// observe
void EpsilonGreedyPolicy::Observe(int a, real r)
{
    estimator->Observe(a, r);
}

/// Destroy
EpsilonGreedyPolicy:: ~EpsilonGreedyPolicy() {}
    
/// Select action: randomly with probability epsilon, otherwise greedily
int EpsilonGreedyPolicy::SelectAction()
{
    if (urandom()<epsilon) {
        return rand()%n_actions;
    }
    return estimator->GetMax();
}


/// Independent population
PopulationOptimalPolicy::PopulationOptimalPolicy(int n_actions, ActionValueEstimate* estimator, real gamma, int n_samples)
{
    this->n_actions = n_actions;
    this->estimator = estimator;
    this->gamma = gamma;
    this->n_samples = n_samples;
}
/// Reset
void PopulationOptimalPolicy::Reset()
{
    estimator->Reset();
}
void PopulationOptimalPolicy::Observe(int a, real r)
{
    estimator->Observe(a, r);
}
PopulationOptimalPolicy::~PopulationOptimalPolicy()
{
}
int PopulationOptimalPolicy::SelectAction()
{
    int j = estimator->GetMax();
    real q_j = estimator->GetMean(j);
    int a = j;
    real U_max = q_j; // Why was it -1 before???
    //real factor = gamma / (1.0f - gamma); // oops, wrong
    real factor = 1.0f / (1.0f - gamma);
    for (int i=0; i<n_actions; i++) {
      if (i==j) continue;
        for (int k=0; k<n_samples; k++) {
            real delta = estimator->Sample(i) - q_j;
            real p = estimator->GetProbability(i, j, delta);
            real U = factor*delta*p - (1-p)*q_j;
            if (U>0) {
                if (U>U_max) {
                    a = i;
                    U_max = i;
                }
            }
        }
    }
    return a;
}

//------- Population Sample Policy ----------//
PopulationSamplePolicy::PopulationSamplePolicy(int n_actions, ActionValueEstimate* estimator, real gamma)
{
    this->n_actions = n_actions;
    this->estimator = estimator;
    this->gamma = gamma;
}
/// Reset
void PopulationSamplePolicy::Reset()
{
    estimator->Reset();
}
void PopulationSamplePolicy::Observe(int a, real r)
{
    estimator->Observe(a, r);
}
PopulationSamplePolicy::~PopulationSamplePolicy() {}
int PopulationSamplePolicy::SelectAction()
{
    return estimator->Sample();
}



/// Create a new E3 policy.
NaiveE3Policy::NaiveE3Policy(int n_actions, real epsilon, real gamma)
{
    this->n_actions = n_actions;
    this->epsilon = epsilon;
    this->gamma = gamma;
    T = 1.0 / (1.0 - gamma);
    estimator = new ActionValueE3Estimate(n_actions, epsilon, T);
}

/// Reset
void NaiveE3Policy::Reset() {
    estimator->Reset();
}
void NaiveE3Policy::Observe(int a, real r)
{
    estimator->Observe(a, r);
}
/// Destroy
NaiveE3Policy:: ~NaiveE3Policy()
{
    delete estimator;
}
    
/// Select action: randomly with probability epsilon, otherwise greedily
int NaiveE3Policy::SelectAction()
{
    return estimator->Sample();
}




/// Make new population VPI policy
PopulationVPIPolicy::PopulationVPIPolicy(int n_actions, PopulationEstimate* estimator, real gamma)
    {
        this->n_actions = n_actions;
        this->estimator = estimator;
        this->gamma = gamma;
    }
    /// Reset
    void PopulationVPIPolicy::Reset() {
        estimator->Reset();
    }
    void PopulationVPIPolicy::Observe(int a, real r)
    {
        estimator->Observe(a, r);
    }
    PopulationVPIPolicy::~PopulationVPIPolicy() {}
    // Each sample is a complete set of samples, not just the one
    // of the best MDP thingy.. so..
    int PopulationVPIPolicy::SelectAction()
    {
        real B = 1.0f / (1.0f - gamma); ///< the exploration factor
        int j = estimator->GetMax();
        int j2 = estimator->GetSecondMax();
        real q_j = B * estimator->GetMean(j);
        real q_j2 = B * estimator->GetMean(j2);
        int a = j;  ///< the selected action
        real U_max = -1;

        // Estimate the utility of all actions
        for (int i=0; i<n_actions; i++) {
            real V = B * estimator->GetMean(i);
            int N = estimator->q[i].N;
            real Gain = 0.0f;
            
            // Ingegrate over gain
            for (int n=0; n<N; n++) {
#if 0
                // Select and MDP and calculate max policy from that
                int sample_argmax = estimator->Sample();
                real sample_max = estimator->s[sample_argmax];
                real q_i = estimator->s[i];
#else
                // Calculate max policy from current value..
                real sample_max = estimator->GetMean(j);
                real q_i = estimator->Sample(i);
#endif
                real q_i_max;
                if (q_i > sample_max) {
                    q_i_max = B*q_i;
                } else {
                    q_i_max = q_i + gamma * B * sample_max;
                }
                if (i==j && q_j2 > q_i_max) {
                    Gain += q_j2 - q_i_max;
                } else if (i!=j && q_i_max > q_j) {
                    Gain += q_i_max - q_j;
                }
            }
            Gain /= (real) N;
            real U = V + Gain;
            if (U > U_max || i == 0) {
                U_max = U;
                a = i;
            }
        }
        return a;
    }



/// A new Particle Filter VPI policy
PFVPIPolicy::PFVPIPolicy(int n_actions, PFActionValueEstimate* estimator, real gamma, int n_samples)
{
    this->n_actions = n_actions;
    this->estimator = estimator;
    this->gamma = gamma;
    this->n_samples = n_samples;
}

/// Reset
void PFVPIPolicy::Reset() {
    estimator->Reset();
}

void PFVPIPolicy::Observe(int a, real r)
{
    estimator->Observe(a, r);
}

PFVPIPolicy::~PFVPIPolicy()
{
}

int PFVPIPolicy::SelectAction()
{
    int j = estimator->GetMax();
    int j2 = estimator->GetSecondMax();
    real q_j = estimator->GetMean(j);
    real q_j2 = estimator->GetMean(j2);
    int a = j;
    real U_max = -1;
    for (int i=0; i<n_actions; i++) {
        real V = estimator->GetMean(i);
        real Gain = 0.0f;
        int N = estimator->q[i].N;
        real sum = 0.0;
        for (int n=0; n<N; n++) {
            real q_i = estimator->q[i].y[n];
            real p_i = estimator->q[i].w[n];
            sum += p_i;
            if (i==j && q_j2 > q_i) {
                Gain += (q_j2 - q_i)*p_i;
            } else if (i!=j && q_i > q_j) {
                Gain += (q_i - q_j)*p_i;
            }
        }
        real U = V + Gain/sum;
        if (U > U_max || i == 0) {
            U_max = U;
            a = i;
        }
    }
    return a;
}

OptimalInfinitePolicy::OptimalInfinitePolicy(int n_actions)
{
    this->n_actions = n_actions;
    estimate = new CountingBernoulliEstimate(n_actions);
}
void OptimalInfinitePolicy::Reset()
{
    estimate->Reset();
}
void OptimalInfinitePolicy::Observe(int a, real r)
{
    estimate->Observe(a, r);
}
OptimalInfinitePolicy::~OptimalInfinitePolicy()
{
    delete estimate;
}
int OptimalInfinitePolicy::SelectAction()
{
    return estimate->Sample();
}
