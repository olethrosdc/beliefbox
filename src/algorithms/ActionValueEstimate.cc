/* -*- Mode: C++; -*- */
/* VER: $Id: ActionValueEstimate.c,v 1.6 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ActionValueEstimate.h"
#include "MathFunctions.h"
#include "Random.h"

PointEstimate::PointEstimate(int n_actions, real alpha, real min, real range)
{
    this->n_actions = n_actions;
    this->alpha = alpha;
    this->min = min;
    this->range = range;

    q.resize(n_actions);
    Reset();
}

/// Reset - to random or something else
void PointEstimate::Reset()
{
    for (int i=0; i<n_actions; i++) {
        q[i] = min + urandom() * range;
    }
}

/// Destroy
PointEstimate::~PointEstimate() {}

/// Update estimates when observing reward r after action a is taken
void PointEstimate::Observe (int a, real r)
{
    q[a] += alpha * (r - q[a]);
}

/// Get the mean reward of action a
real PointEstimate::GetMean (int a)
{
    return q[a];
}

/// Sample an action i with probability p(q[i] > q[j])
int PointEstimate::Sample()
{
    return GetMax();
}

/// Sample the action with maximum mean
int PointEstimate::GetMax()
{
    int arg_max = 0;
    real max = q[0];
    for (int i=1; i<n_actions; i++) {
        if (q[i] > max) {
            arg_max = i;
            max = q[i];
        }
    }
    return arg_max;
}


int PointEstimate::GetSecondMax()
{
    int arg_max = GetMax();
    int arg_max2 = 0;
    real max;
    if (arg_max==0) {
        max = q[1];
    } else {
        max = q[0];
    }
    for (int i=0; i<n_actions; i++) {
        if (i!=arg_max){
            real m = q[i];
            if (max < m) {
                max = m;
                arg_max2 = i;
            }
        }
    }
    return arg_max2;
}

real PointEstimate::Sample(int a)
{
    return q[a];
}

PopulationEstimate::PopulationEstimate(int n_actions, int n_members, real alpha)
{
    this->n_actions = n_actions;
    q.resize(n_actions);
    s.resize(n_actions);
    for (int i=0; i<n_actions; i++) {
        q[i].SetNumberOfEstimates (n_members);
        q[i].SetLearningRate (alpha);
        q[i].Reset();
    }
}

void PopulationEstimate::Reset()
{
    for (int i=0; i<n_actions; i++) {
        q[i].Reset();
    }
}

/// Destroy
PopulationEstimate::~PopulationEstimate()
{
}

/// Update estimates when observing reward r after action a is taken
void PopulationEstimate::Observe (int a, real r)
{
    q[a].Observe(r);
}

/// Get the mean reward of action a
real PopulationEstimate::GetMean (int a)
{
    return q[a].GetMean();
}

/// Get the mean reward of action a
real PopulationEstimate::GetVar (int a)
{
    return q[a].GetVar();
}

/// Sample an action i with probability p(q[i] > q[j])
int PopulationEstimate::Sample()
{
    for (int i=0; i<n_actions; i++) {
        s[i] = q[i].Sample();
    }
    int arg_max = 0;
    real max = s[0];
    for (int i=1; i<n_actions; i++) {
        if (max < s[i]) {
            max = s[i];
            arg_max = i;
        }
    }
    return arg_max;
}
/// Sample the action with maximum mean
int PopulationEstimate::GetMax()
{
    int arg_max = 0;
    real max = q[0].GetMean();
    for (int i=1; i<n_actions; i++) {
        real m = q[i].GetMean();
        if (max < m) {
            max = m;
            arg_max = i;
        }
    }
    return arg_max;
}
int PopulationEstimate::GetSecondMax()
{
    int arg_max = GetMax();
    int arg_max2 = 0;
    real max;
    if (arg_max==0) {
        max = q[1].GetMean();
    } else {
        max = q[0].GetMean();
    }
    for (int i=0; i<n_actions; i++) {
        if (i!=arg_max){
            real m = q[i].GetMean();
            if (max < m) {
                max = m;
                arg_max2 = i;
            }
        }
    }
    return arg_max2;
}

real PopulationEstimate::Sample(int a)
{
    return q[a].Sample();
}


real PopulationEstimate::GetMemberValue(int a, int i)
{
    return q[a].x[i];
}

/// Get a simple estimate of P(q_i - q_j > delta)
real PopulationEstimate::GetProbability(int i, int j, real delta)
{
    real P = 0.0f;
    int N = q[i].N;
    for (int n=0; n<N; n++) {
        if (q[i].x[n] - q[j].x[n] > delta) {
            P += 1.0f;
        }
    }
    return P / (real) N;
}


AverageEstimate::AverageEstimate ()
{
    init = 1.0f;
    Reset();
}

AverageEstimate::AverageEstimate (real init_, int N_)
{
    init = init_;
    N_init = N_;
    Reset();
}

void AverageEstimate::Reset() {
    Ex = init;
    N = N_init;
}

void AverageEstimate::Observe (real X)
{
    real fN = (real) N;
    Ex = (Ex*fN + X)/(fN + 1.0);
    N++;
}

real AverageEstimate::GetMean()
{
    return Ex;
}

int AverageEstimate::GetN()
{
    return N;
}


/// This is a fully Bayesian estimate, yay
BernoulliEstimate::BernoulliEstimate(int n_actions, int n_samples, real alpha, real beta)
{
    this->n_actions = n_actions;
    this->n_samples = n_samples;
    this->alpha = alpha;
    this->beta = beta;
    Reset();
}

/// Destroy
BernoulliEstimate::~BernoulliEstimate() {}

/// Update estimates when observing reward r after action a is taken
void BernoulliEstimate::Observe (int a, real r)
{
    prior[a].calculatePosterior(r);
}

/// Get the mean reward of action a
real BernoulliEstimate::GetMean (int a)
{
    return prior[a].getMean();
}

/// Sample an action i with probability p(q[i] > q[j])
int BernoulliEstimate::Sample ()
{
    int arg_max = 0;
    real max = prior[0].generate();
    for (int i=1; i<n_actions; i++) {
        real m = prior[i].generate();
        if (max<m) {
            max = m;
            arg_max = i;
        }
    }
    return arg_max;
}

/// Sample the action with maximum mean
int BernoulliEstimate::GetMax ()
{
    int arg_max = 0;
    real max = prior[0].getMean();
    for (int i=1; i<n_actions; i++) {
        real m = prior[i].getMean();
        if (max < m) {
            max = m;
            arg_max = i;
        }
    }
    return arg_max;
}

int BernoulliEstimate::GetSecondMax()
{
    int arg_max = GetMax();
    int arg_max2 = 0;
    real max;
    if (arg_max==0) {
        max = prior[1].getMean();
    } else {
        max = prior[0].getMean();
    }

    for (int i=0; i<n_actions; i++) {
        if (i!=arg_max){
            real m = prior[i].getMean();
            if (max < m) {
                max = m;
                arg_max2 = i;
            }
        }
    }
    return arg_max2;
}

/// Sample from the mean of an action a
real BernoulliEstimate::Sample (int a)
{
    return prior[a].generate();
}
/// Reset estimates
void BernoulliEstimate::Reset ()
{
    prior.resize(n_actions);
    for (int i=0; i<n_actions; i++) {
        prior[i].alpha = alpha;
        prior[i].beta = beta;
    }
}

/// Get a simple estimate of P(q_i - q_j > delta)
real BernoulliEstimate::GetProbability(int i, int j, real delta) 
{
    int N = 0;
    for (int k=0; k<n_samples; k++) {
        real q_i = prior[i].generate();
        real q_j = prior[j].generate();
        if (q_i - q_j > delta) {
            N++;
        }
    }
    return (real) N / (real) n_samples;
}


ActionValueE3Estimate::ActionValueE3Estimate (int n_actions, real epsilon, real T)
{
    this->epsilon = epsilon;
    this->T = T;
    this->n_actions = n_actions;
    q.resize(n_actions);
    Reset();
}
/// Destroy
ActionValueE3Estimate::~ActionValueE3Estimate() {}

/// Update estimates when observing reward r after action a is taken
void ActionValueE3Estimate::Observe (int a, real r)
{
    q[a].Observe(r);
    n_visits++;
}

/// Get the mean reward of action a
real ActionValueE3Estimate::GetMean (int a)
{
    return q[a].GetMean();
}

/// Sample an action i with balance wandering
int ActionValueE3Estimate::Sample () 
{
    real m_known = pow(T*((real) n_actions)/epsilon, 4.0);
    //real m_known = T*((real) n_actions)/epsilon;
    if (n_visits<(int) m_known) {
        // Do balanced wandering
        int arg_min = 0;
        int min = q[0].GetN();
        for (int i=1; i<n_actions; i++) {
            int n = q[i].GetN();
            if (min > n) {
                min = n;
                arg_min = i;
            }
        }
        return arg_min;
    }

    // Be greedy
    return GetMax();
}

/// Sample the action with maximum mean
int ActionValueE3Estimate::GetMax()
{
    int arg_max = 0;
    real max = q[0].GetMean();
    for (int i=1; i<n_actions; i++) {
        real m = q[i].GetMean();
        if (max < m) {
            max = m;
            arg_max = i;
        }
    }
    return arg_max;
}

/// Sample the action with second maximum mean
int ActionValueE3Estimate::GetSecondMax()
{
    int arg_max = GetMax();
    int arg_max2 = 0;
    real max;
    if (arg_max==0) {
        max = q[1].GetMean();
    } else {
        max = q[0].GetMean();
    }
    for (int i=0; i<n_actions; i++) {
        if (i!=arg_max){
            real m = q[i].GetMean();
            if (max < m) {
                max = m;
                arg_max2 = i;
            }
        }
    }
    return arg_max2;
}

/// Sample from an action i  - currently returns mean;
real ActionValueE3Estimate::Sample (int a)
{
    return q[a].GetMean();
}
/// Reset estimates
void ActionValueE3Estimate::Reset ()
{
    n_visits = 0;
    for (int i=0; i<n_actions; i++) {
        q[i].Reset();
    }
}

CountingBernoulliEstimate::CountingBernoulliEstimate(int n_actions)
{
    this->n_actions = n_actions;
    obs = new BernoulliArmObservations [n_actions];
    Reset();
}

CountingBernoulliEstimate::~CountingBernoulliEstimate()
{
    delete [] obs;
}

/// Update estimates when observing reward r after action a is taken
void CountingBernoulliEstimate::Observe (int a, real r)
{
    if (r>0.99) {
        obs[a].n_succ++;
    } else {
        obs[a].n_fail++;
    }
    //printf ("O: %d %f, ", a, r);
}

/// Get the mean reward of action a
real CountingBernoulliEstimate::GetMean (int a)
{
    return (float) obs[a].n_succ / (float) (obs[a].n_succ + obs[a].n_fail);
}
/// Sample an action i with probability p(q[i] > q[j])
int CountingBernoulliEstimate::Sample ()
{
  int a = ArgMax (n_actions, obs);
  //  for (int i=0; i<n_actions; i++) {
  //  printf ("(%d: %d %d)", i, obs[i].n_succ, obs[i].n_fail);
  //}
  //printf (" %d\n", a);  
  return a;
}
/// Sample the action with maximum mean
int CountingBernoulliEstimate::GetMax ()
{
    return Sample();
}
/// Sample the second action
int CountingBernoulliEstimate::GetSecondMax()
{
    return Sample();
}

/// Sample an action i with probability p(q[i] > q[j])
real CountingBernoulliEstimate::Sample (int a)
{
    return GetMean(a);
}
/// Reset estimates
void CountingBernoulliEstimate::Reset ()
{
    for (int i=0; i<n_actions; i++) {
        obs[i].n_succ = 0;
        obs[i].n_fail = 0;
    }
}

real CountingBernoulliEstimate::GetProbability(int i, int j, real delta) 
{
    return 0.0f;
}

