/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_BANDIT_POLICY_H
#define DISCRETE_BANDIT_POLICY_H

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "SampleEstimator.h"
#include "ActionValueEstimate.h"
#include "PFActionValueEstimate.h"

/** A general policy
 */
class DiscreteBanditPolicy 
{
public:
    /// Destructor
    virtual ~DiscreteBanditPolicy() {}
    /// Select an action according to this policy
    virtual int SelectAction() = 0;
    virtual void Reset() = 0;
    virtual void Observe(int a, real r) = 0;
};



/** A greedy policy
 */
class EpsilonGreedyPolicy : public DiscreteBanditPolicy
{
public:
    int n_actions;
    real epsilon;
    ActionValueEstimate* estimator;
    EpsilonGreedyPolicy(int n_actions, real epsilon, ActionValueEstimate* estimator);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~EpsilonGreedyPolicy();
    virtual int SelectAction();
};


/** An optimal policy.

This policy is optimistically optimal.
*/
class PopulationOptimalPolicy : public DiscreteBanditPolicy
{
 public:
    int n_actions;
    ActionValueEstimate* estimator;
    real gamma;
    int n_samples;
    PopulationOptimalPolicy(int n_actions, ActionValueEstimate* estimator, real gamma, int n_samples);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~PopulationOptimalPolicy();
    virtual int SelectAction();
};



/** A sampling optimal policy.

This policy is asymptoticall optimal.
*/
class PopulationSamplePolicy : public DiscreteBanditPolicy
{
 public:
    int n_actions;
    ActionValueEstimate* estimator;
    real gamma;
    PopulationSamplePolicy(int n_actions, ActionValueEstimate* estimator, real gamma);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~PopulationSamplePolicy();
    virtual int SelectAction();
};




/** A naive version of the E3 algorithm
    
 We have an MDP \f$M\f$ with states \f$\{1, .., N\}\f$ and action
 \f$\{a_i\}_{i=1}^k\f$.

- \f$P_M^a(ij) \geq 0\f$, transitions from \f$i\f$ to \f$j\f$.
- \f$R_(i), R_{max} \geq R_M(i) \geq 0, Var_M(i) \leq Var_{max}\f$

I disagree with the discussion that there is no 'unambiguously' better policy.

- T-path in M:  a sequence p of T+1 states,
\f$p = \{i_1, i_2, \cdots, i_{T+1}\}\f$

- Return along p \f$U_M(p) = (1/T)(R_{i_1}+ \ldots +R_{i_T})\f$

- T-step return from state i \f$U_M^\pi(p) = \sum_p P_M^\pi[p]U_m(p)\f$,
where the sum is over all T-paths p in M that start at i.

- \f$U_M^\pi(i) = \lim_{T \to \infty} U_M^\pi(i,T)\f$.  Since we
  are in the unichain, we are independent of i, so we just write
  \f$U_M^\pi\f$ (but this is not the case for discounting.

- Optimal T-step return \f$U_M^*(i,T) = \max_\pi \{U_M^\pi(i,T)\}\f$.
  Also \f$U_M^*(i) = \lim_{T \to \infty} U_M^*(i,T)\f$.

- The maximum possible T-step return is \f$R_{max}\f$.

- The \f$\epsilon\f$-return mixing time of \f$\pi\f$ is the smallest \f$T\f$ such that for all \f$T' \geq T, |U_M^\pi(i,T') - U_M^\pi| \leq \epsilon\f$ for all \f$i\f$.



 Theorem 1: There exists an algortihm \f$A\f$, taking inputs
\f$\epsilon, \delta, N, T, opt(\Pi_M^{T,\epsilon})\f$, such that if
the total number of actions and computation time taken by A exceeds a
polynomial in \f$1/\epsilon, 1/\delta, N, T, R_{max}\f$, then with
probability at least \f$1-\delta\f$, the total undescounted return of \f$A\f$ will exceed \f$opt(\Pi_M^{T,\epsilon}) - \epsilon\f$.

 */
class NaiveE3Policy: public DiscreteBanditPolicy
{
public:
    int n_actions;
    real epsilon;
    real gamma;
    real T;
    ActionValueE3Estimate* estimator;

    /// Create a new e-greedy policy
    NaiveE3Policy(int n_actions, real epsilon, real gamma);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~NaiveE3Policy();
    virtual int SelectAction();
};



/** A VPI policy for population estimates.

This policy is optimal according to the VPI criterion.
It is specialised to a population estimate.
*/
class PopulationVPIPolicy : public DiscreteBanditPolicy
{
 public:
    int n_actions;
    PopulationEstimate* estimator;
    real gamma;

    PopulationVPIPolicy(int n_actions, PopulationEstimate* estimator, real gamma);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~PopulationVPIPolicy();
    virtual int SelectAction();
};

/** A VPI policy.

This policy is optimal according to the VPI criterion.
It is not specialised for a particular estimate.
*/
class VPIPolicy : public DiscreteBanditPolicy
{
 public:
    int n_actions;
    ActionValueEstimate* estimator;
    real gamma;
    int n_samples;
    VPIPolicy(int n_actions, ActionValueEstimate* estimator, real gamma, int n_samples);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~VPIPolicy();
    virtual int SelectAction();
};


/** A VPI policy for particle filter estimates.

This policy is optimal according to the VPI criterion.
It is specialised for a particle filter estimate.
*/
class PFVPIPolicy : public DiscreteBanditPolicy
{
 public:
    int n_actions;
    PFActionValueEstimate* estimator;
    real gamma;
    int n_samples;
    PFVPIPolicy(int n_actions, PFActionValueEstimate* estimator, real gamma, int n_samples);
    /// Reset
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~PFVPIPolicy();
    virtual int SelectAction();
};


class OptimalInfinitePolicy : public DiscreteBanditPolicy
{
public:
    int n_actions;
    CountingBernoulliEstimate* estimate;
    OptimalInfinitePolicy(int n_actions);
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual ~OptimalInfinitePolicy();
    virtual int SelectAction();
};



#endif
