// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: ValueIteration.c,v 1.5 2006/11/08 17:20:17 cdimitrakakis Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MultiMDPValueIteration.h"
#include "real.h"
#include "MathFunctions.h"
#include "Vector.h"
#include <cmath>
#include <cassert>

/// Setup a value iteration procedure for a list of MDPs with weights w and discount factor gamma
MultiMDPValueIteration::MultiMDPValueIteration(const Vector& w_,
                                               const std::vector<const DiscreteMDP*> &mdp_list_,
                                               real gamma_) :
    w(w_),
    mdp_list(mdp_list_),
    gamma(gamma_),
    n_mdps(mdp_list.size())
{
    assert (mdp_list.size());
    assert (gamma>=0 && gamma <=1);
    assert ((int) mdp_list.size() == w.Size());
    n_actions = mdp_list[0]->GetNActions();
    n_states = mdp_list[0]->GetNStates();
    for (int i=1; i<n_mdps; ++i) {
        if (n_actions != mdp_list[i]->GetNActions()) {
            throw std::runtime_error("Number of actions in MDPs does not agree\n");
        }
        if (n_states != mdp_list[i]->GetNStates()) {
            throw std::runtime_error("Number of states in MDPs does not agree\n");
        }
    }
    Reset();
}

/// Reset
void MultiMDPValueIteration::Reset()
{
    V.Resize(n_states);
    dV.Resize(n_states);
    pV.Resize(n_states);
    Q.Resize(n_states, n_actions);
    dQ.Resize(n_states, n_actions);
    pQ.Resize(n_states, n_actions);
    for (int s=0; s<n_states; s++) {
        V(s) = 0.0;
        dV(s) = 0.0;
        pV(s) = 0.0;
        for (int a=0; a<n_actions; a++) {
            Q(s,a) = 0.0;
            dQ(s,a) = 0.0;
            pQ(s,a) = 0.0;
        }
    }
    V_mu.resize(n_mdps);
    for (int i=0; i<n_mdps; ++i) {
        V_mu[i].Resize(n_states);
        V_mu[i].Clear();
    }
}


/// Empty destructor
MultiMDPValueIteration::~MultiMDPValueIteration()
{
}


/** Return the state-action value for the next stage.

    The function returns
    \f[
    Q_{\mu,t}(s,a) = r_\mu(s,a) + \gamma \sum_{s'} V_{\mu, t+1}(s') {\cal T}_\mu(s' | s, a)
    \f]
 */
real MultiMDPValueIteration::ComputeStateActionValueForSingleMDP(int mu, int s, int a)
{
    Vector& V_i = V_mu[mu];
    const DiscreteMDP* mdp = mdp_list[mu];
    real Q_mu_sa = mdp->getExpectedReward(s,a);
    const DiscreteStateSet& next = mdp->getNextStates(s, a);
    for (DiscreteStateSet::iterator i=next.begin();
         i!=next.end();
         ++i) {
        int s2 = *i;
        real P = mdp->getTransitionProbability(s, a, s2);
        real R = gamma * V_i(s2);
        Q_mu_sa += P * R;
    } 
    return Q_mu_sa;
}

/*** Compute the next state values for the MDPs.
     
     This is supposed to be used internally by the value iteration.
     It calculates 
     \f[
     Q_{\xi, t}(s,a) = \sum_\mu \xi(\mu) Q_{\mu, t}(s,a)
     \f]
     
     with

     \f[
      Q_{\mu, t}(s,a) = \sum_{s'} {\cal T}_\mu (s'|s,a) V_{\mu, t+1}(s')
     \f]

 */
real MultiMDPValueIteration::ComputeActionValueForMDPs(int s, int a)
{
    real Q_sa = 0;
    for (int i=0; i<n_mdps; ++i) {
        real Q_mu_sa = ComputeStateActionValueForSingleMDP(i, s, a);
        Q_sa += w(i) * Q_mu_sa;
    } 
    return Q_sa;
}

/*** Compute the current value.
     
     First, we require the calculation of the optimal policy for all
     states for the current step.  This is given by:
     
     \[f
     a^*_{\xi, t}(s) = \arg\max_a Q_{\xi,t}(s,a)
     \f]
 */
void MultiMDPValueIteration::ComputeStateValues(real threshold, int max_iter)
{
        //Vector pV(V.size());
        //Vector dV(V.size());
    pV = V;
    do {
        Delta = 0.0;
        std::vector<int> a_max(n_states);
        for (int s=0; s<n_states; s++) {
            real Q_a_max = -RAND_MAX;
            a_max[s] = 0;
            for (int a=0; a<n_actions; a++) {
                real Q_sa = ComputeActionValueForMDPs(s, a);
                if (a==0 || Q_a_max < Q_sa) {
                    a_max[s] = a;
                    Q_a_max = Q_sa;
                }
            }
            V(s) = Q_a_max;
        }
        for (int i=0; i<n_mdps; ++i) {
            Vector tmpV = V_mu[i];
            for (int s=0; s<n_states; s++) {
                tmpV(s) = ComputeStateActionValueForSingleMDP(i, s, a_max[s]);
            }
            V_mu[i] = tmpV;
        }
        if (max_iter > 0) {
            max_iter--;
        }
        Delta = abs(V - pV).Sum();
		printf("D:%f i:%d\n", Delta, max_iter);
    } while(Delta >= threshold && max_iter != 0);
	
}




/** ComputeStateActionValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/

void MultiMDPValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
         //Vector pV(V.size());
        //Vector dV(V.size());
    int n_iter = 0;
    do {
        pV = V;
        Delta = 0.0;

        // Find the  best average action at the current stage.
        std::vector<int> a_max(n_states); 
        for (int s=0; s<n_states; s++) {
            real Q_a_max = -RAND_MAX;
            a_max[s] = 0;
            for (int a=0; a<n_actions; a++) {
                real Q_sa = ComputeActionValueForMDPs(s, a);
                if (a==0 || Q_a_max < Q_sa) {
                    a_max[s] = a;
                    Q_a_max = Q_sa;
                }
                Q(s,a) = Q_sa;
            }
            V(s) = Q_a_max;
        }

        // Calculate the value of each MDP for the current best average action
        for (int i=0; i<n_mdps; ++i) {
            Vector tmpV = V_mu[i];
            for (int s=0; s<n_states; s++) {
                tmpV(s) = ComputeStateActionValueForSingleMDP(i, s, a_max[s]);
            }
            V_mu[i] = tmpV;
        }
        if (max_iter > 0) {
            max_iter--;
        }
        Delta = abs(V - pV).Sum();
        n_iter++;
    } while(Delta >= threshold && max_iter != 0);
    //printf("Exiting at d:%f, n:%d\n", Delta, n_iter);	
}

/// Get the policy that is optimal for the mixed MDP.
FixedDiscretePolicy* MultiMDPValueIteration::getPolicy()
{
    FixedDiscretePolicy* policy = new FixedDiscretePolicy(n_states, n_actions);
    for (int s=0; s<n_states; s++) {
        real max_Qa = getValue(s, 0);
        int argmax_Qa = 0;
        for (int a=1; a<n_actions; a++) {
            real Qa = getValue(s, a);
            if (Qa > max_Qa) {
                max_Qa = Qa;
                argmax_Qa = a;
            }
        }
        Vector* p = policy->getActionProbabilitiesPtr(s);
        for (int a=0; a<n_actions; a++) { 
            (*p)[a] = 0.0;
        }
        (*p)[argmax_Qa] = 1.0;
    }
    return policy;
}


