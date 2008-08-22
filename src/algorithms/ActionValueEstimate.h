/* -*- Mode: C++; -*- */
/* VER: $Id: ActionValueEstimate.h,v 1.4 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ACTION_VALUE_ESTIMATE
#define ACTION_VALUE_ESTIMATE

#include <cmath>
#include <vector>

#include "real.h"
#include "SampleEstimator.h"
#include "ParticleFilter.h"
#include "BetaDistribution.h"

/** An Estimate of actions.
 */
class ActionValueEstimate
{
 public:
    /// Destroy
    virtual ~ActionValueEstimate() {}
    /// Update estimates when observing reward r after action a is taken
    virtual void Observe (int a, real r) = 0;
    /// Get the mean reward of action a
    virtual real GetMean (int a) = 0;
    /// Sample an action i with probability p(q[i] > q[j])
    virtual int Sample () = 0;
    /// Sample the action with maximum mean
    virtual int GetMax () = 0;
    /// Sample the second action
    virtual int GetSecondMax() = 0;
    /// Sample an action i with probability p(q[i] > q[j])
    virtual real Sample (int a) = 0;
    /// Reset estimates
    virtual void Reset () = 0;
    virtual real GetProbability(int i, int j, real delta) 
    {
        return 0.0f;
    }
};


/** A point estimate of actions
 */
class PointEstimate : public ActionValueEstimate
{
 public:
    std::vector<real> q;
    int n_actions;
    real alpha;
    real min;
    real range;
    PointEstimate(int n_actions, real alpha, real min=0.0, real range = 1.0);
    /// Reset - to random or something else
    virtual void Reset();
    virtual ~PointEstimate();
    /// Update estimates when observing reward r after action a is taken
    virtual void Observe (int a, real r);
    /// Get the mean reward of action a
    virtual real GetMean (int a);
    /// Sample an action i with probability p(q[i] > q[j])
    virtual int Sample();
    /// Sample the action with maximum mean
    virtual int GetMax();
    virtual int GetSecondMax();
    virtual real Sample(int a);
};

/** A population estimate of actions
 */
class PopulationEstimate : public ActionValueEstimate
{
public:
    std::vector<SampleEstimator> q;
    std::vector<real> s;
    int n_actions;
    PopulationEstimate(int n_actions, int n_members, real alpha);
    virtual void Reset();
    /// Destroy
    virtual ~PopulationEstimate();
    /// Update estimates when observing reward r after action a is taken
    virtual void Observe (int a, real r);
    /// Get the mean reward of action a
    virtual real GetMean (int a);
    /// Get the mean reward of action a
    virtual real GetVar (int a);
    /// Sample an action i with probability p(q[i] > q[j])
    virtual int Sample();
    /// Sample the action with maximum mean
    virtual int GetMax();
    virtual int GetSecondMax();
    virtual real Sample(int a);
    virtual real GetMemberValue(int a, int i);
    /// Get a simple estimate of P(q_i - q_j > delta)
    virtual real GetProbability(int i, int j, real delta);
};

class AverageEstimate
{
public:
    real Ex;
    int N;
    real init;
    int N_init;
    AverageEstimate ();
    AverageEstimate (real init_, int N_ = 0);
    void Reset();
    void Observe (real X);
    real GetMean();
    int GetN();
};


/// This is a fully Bayesian estimate, yay
class BernoulliEstimate : public ActionValueEstimate
{
 public:
    int n_actions; ///< number of actions to estimate
    int n_samples; ///< number of samples to take for estimation
    real alpha;
    real beta;
    std::vector<BetaDistribution> prior;
    BernoulliEstimate(int n_actions, int n_samples, real alpha=1.0, real beta=1.0);
    /// Destroy
    virtual ~BernoulliEstimate();
    /// Update estimates when observing reward r after action a is taken
    virtual void Observe (int a, real r);
    /// Get the mean reward of action a
    virtual real GetMean (int a);
    /// Sample an action i with probability p(q[i] > q[j])
    virtual int Sample ();
    /// Sample the action with maximum mean
    virtual int GetMax ();
    virtual int GetSecondMax();
    /// Sample an action i with probability p(q[i] > q[j])
    virtual real Sample (int a);
    /// Reset estimates
    virtual void Reset ();
    virtual real GetProbability(int i, int j, real delta);
};


/** Action value estimate for E3
 */
class ActionValueE3Estimate : public ActionValueEstimate
{
 public:
    real epsilon;
    real T;
    int n_actions;
    int n_visits;
    std::vector<AverageEstimate> q;
    ActionValueE3Estimate (int n_actions, real epsilon, real T);
    /// Destroy
    virtual ~ActionValueE3Estimate();
    /// Update estimates when observing reward r after action a is taken
    void Observe (int a, real r);
    /// Get the mean reward of action a
    virtual real GetMean (int a);
    /// Sample an action i with balance wandering
    virtual int Sample ();
    /// Sample the action with maximum mean
    virtual int GetMax();
    virtual int GetSecondMax();
    /// Sample from an action i  - currently returns mean;
    virtual real Sample (int a);
    /// Reset estimates
    virtual void Reset ();
};


class BernoulliArmObservations
{
public:
    int n_succ;
    int n_fail;
    BernoulliArmObservations()
    {
        n_succ = 0;
        n_fail = 0;
    }

    /// Compare two times
    bool operator< (const BernoulliArmObservations& rhs) const
    {
        if (n_fail > rhs.n_fail) {
            return true;
        } else if (n_fail < rhs.n_fail) {
            return false;
        } else if (n_succ < rhs.n_succ) {
            return true;
        }
        return false;
    }
    /// Compare two times
    bool operator> (const BernoulliArmObservations& rhs) const
    {
        if (n_fail < rhs.n_fail) {
            return true;
        } else if (n_fail > rhs.n_fail) {
            return false;
        } else if (n_succ > rhs.n_succ) {
            return true;
        }
        return false;
    }

};


/** An Estimate of actions.
 */
class CountingBernoulliEstimate : public ActionValueEstimate
{
public:
    BernoulliArmObservations* obs;
    int n_actions;
    CountingBernoulliEstimate(int n_actions);
    /// Destroy
    virtual ~CountingBernoulliEstimate();
    /// Update estimates when observing reward r after action a is taken
    virtual void Observe (int a, real r);
    /// Get the mean reward of action a
    virtual real GetMean (int a);
    /// Sample an action i with probability p(q[i] > q[j])
    virtual int Sample ();
    /// Sample the action with maximum mean
    virtual int GetMax ();
    /// Sample the second action
    virtual int GetSecondMax();
    /// Sample an action i with probability p(q[i] > q[j])
    virtual real Sample (int a);
    /// Reset estimates
    virtual void Reset ();
    virtual real GetProbability(int i, int j, real delta);
};

#endif
