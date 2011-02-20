// -*- Mode: C++ -*-

#ifndef MWAL_H
#define MWAL_H

#include "Vector.h"
#include "Demonstrations.h"
#include "DiscreteMDP.h"
#include "DiscretePolicy.h"


/** Multiplicative weights for apprenticeship learning (MWAL) algorithm.

 Inputs: MDP M, Feature expectations \f$mu_E\f$.

 The feature expectations \f$mu_E\f$ is just an estimate of the
 discounted state distribution.

 The assumption is that there is some \f$w^*\f$ such that \f$R^*(s) = w^* \cdot
 \phi(s)\f$ for some \f$w^*\f$, with \f$\|w\|_1 = , w \geq 0\f$. We wish to find
 \f[
 v^* = max_\psi min_w [w \cdot \mu(\psi) - w \mu_E].
 \f]
 We will find a mixed polciy \f$\psi^*\f$ approximating \f$v^*\f$. Note that
 \f$V(\psi) = w^* \mu(\psi)\f$.
 Consequenty, \f$\psi^*\f$ maximises \f$V(\psi) -
 V(\pi_E)\f$, with respect to the worst-case \f$w\f$.

 The form of the above is the same as a two-player zero-sum game. with
 game matrix:
 \f[
 G(i, j) = \mu^j(i) - \mu_E(i),
 \f]
 where \f$\mu^j\f$ is the vector of feature expecations for
 the \f$j\f$-th policy.

 If there are \f$S\f$ states and \f$A\f$ actions, that means there are
 \f$A^S\f$ deterministic policies, but we do not need to enumerate
 them!

 policy is a deterministic policy, expressed as a mapping S -> A
 mean_policy is a stochastic policy, expressed as a SxA matrix.
 W is the opponent's best response.
*/
class MWAL
{
public:
    int n_states;
    int n_actions;
    real gamma;
    Vector mu_E; ///< discounted feature counts
    FixedDiscretePolicy mean_policy; ///< The policy
    MWAL(int n_states_, int n_actions_, real discount) 
        : n_states(n_states_), n_actions(n_actions_), gamma(discount),
          mu_E(n_states), mean_policy(n_states, n_actions)
    {}

    void CalculateFeatureCounts(Demonstrations<int, int>& D);
    Vector CalculateFeatureExpectation(DiscreteMDP& mdp, 
                                       FixedDiscretePolicy& policy, 
                                       real gamma, real epsilon);
    void Compute(DiscreteMDP& mdp_copy, real gamma, real epsilon, int T = -1);
};

#endif
