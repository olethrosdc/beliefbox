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
    n_actions = mdp_list[0]->getNActions();
    n_states = mdp_list[0]->getNStates();
    for (int i=1; i<n_mdps; ++i) {
        if (n_actions != mdp_list[i]->getNActions()) {
            throw std::runtime_error("Number of actions in MDPs does not agree\n");
        }
        if (n_states != mdp_list[i]->getNStates()) {
            throw std::runtime_error("Number of states in MDPs does not agree\n");
        }
    }
    Reset();
}

/// Reset
void MultiMDPValueIteration::Reset()
{
    V_xi.Resize(n_states);
    dV_xi.Resize(n_states);
    pV_xi.Resize(n_states);
    Q_xi.Resize(n_states, n_actions);
    dQ_xi.Resize(n_states, n_actions);
    pQ_xi.Resize(n_states, n_actions);
    for (int s=0; s<n_states; s++) {
        V_xi(s) = 0.0;
        dV_xi(s) = 0.0;
        pV_xi(s) = 0.0;
        for (int a=0; a<n_actions; a++) {
            Q_xi(s,a) = 0.0;
            dQ_xi(s,a) = 0.0;
            pQ_xi(s,a) = 0.0;
        }
    }
    V.resize(n_mdps);
    Q.resize(n_mdps);
    for (int i=0; i<n_mdps; ++i) {
        V[i].Resize(n_states);
        V[i].Clear();
        Q[i].Resize(n_states, n_actions);
        Q[i].Clear();
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
    Vector& V_i = V[mu];
    const DiscreteMDP* mdp = mdp_list[mu];
    real Q_mu_sa = 0.0; 
    const DiscreteStateSet& next = mdp->getNextStates(s, a);
    for (DiscreteStateSet::iterator i=next.begin();
         i!=next.end();
         ++i) {
        int s2 = *i;
        real P = mdp->getTransitionProbability(s, a, s2);
        real V2 = V_i(s2);
        Q_mu_sa += P * V2;
    } 
    real r_sa = mdp->getExpectedReward(s,a);
    Q_mu_sa = r_sa + gamma * Q_mu_sa;
    Q[mu](s,a) = Q_mu_sa;
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
    for (int mu=0; mu<n_mdps; ++mu) {
        real Q_mu_sa = ComputeStateActionValueForSingleMDP(mu, s, a);
        Q_sa += w(mu) * Q_mu_sa;
    } 
    return Q_sa;
}

/*** Compute the current value.

     This is an alternative implementation with no helper function calls.
     Should work the same.

     First, we require the calculation of the optimal policy for all
     states for the current step.  This is given by:
     
     \[f
     a^*_{\xi, t}(s) = \arg\max_a Q_{\xi,t}(s,a)
     \f]
 */
void MultiMDPValueIteration::ComputeStateValues(real threshold, int max_iter)
{
    pV_xi = V_xi;
    int n_iter = 0;

    //logmsg ("Runnign ComputeStateValues with epsilon: %f, iter: %d, gamma: %f", threshold, max_iter, gamma);
    do {
        // Calculate Q_{mu,t}(s,a) from V_mu(s)
        for (int mu=0; mu<n_mdps; ++mu) {
            const DiscreteMDP* mdp = mdp_list[mu];
            for (int s=0; s<n_states; ++s) {
                for (int a=0; a<n_actions; ++a) {
                    real R_sa = mdp->getExpectedReward(s, a);
                    real Q_msa = 0.0;
                    for (int s2=0; s2<n_states; ++s2) {
                        real P = mdp->getTransitionProbability(s, a, s2);
                        real V2 = V[mu](s2);
                        Q_msa += P * V2;
                    }
                    Q[mu](s, a) = R_sa + gamma * Q_msa;
                    //printf ("Q_%d(%d, %d) = %f = %f + %f * %f\n",
                    //mu, s, a, Q[mu](s,a),
                    //R_sa, gamma, Q_msa);
                }
            }
            
        }

        // calculate Q_{xi, t}(s,a)
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                Q_xi(s,a) = 0.0;
                for (int mu=0; mu<n_mdps; ++mu) {
                    Q_xi(s,a) += w(mu) * Q[mu](s,a);
                }
                //printf ("Q_xi(%d, %d) = %f\n", s, a, Q_xi(s,a));
            }
        }
        

        // Calculate a_t^*(s), V_{xi, t}(s)
        std::vector<int> a_max(n_states);         
        for (int s=0; s<n_states; ++s) {
            a_max[s] = 0;
            real Q_max = Q_xi(s,0);
            for (int a=1; a<n_actions; ++a) {
                if (Q_max < Q_xi(s, a)) {
                    Q_max = Q_xi(s, a);
                    a_max[s] = a;
                }
            }
            V_xi(s) = Q_max;
            //printf ("V_xi(%d) = %f\n", s, V_xi(s));
        }

        // Calculate V_{mu, t}(s)
        for (int mu=0; mu<n_mdps; ++mu) {
            for (int s=0; s<n_states; ++s) {
                V[mu](s) = Q[mu](s, a_max[s]);
            }
        }


        Delta = 0.0;
        if (max_iter > 0) {
            max_iter--;
        }
        //V.print(stdout);
        //pV.print(stdout);
        Delta = abs(V_xi - pV_xi).Sum();
        pV_xi = V_xi;
		//printf("%f # delta\n", Delta);
        n_iter++;
    } while(Delta >= threshold && max_iter != 0);
    //logmsg("Exiting at delta :%f, iter :%d\n", Delta, n_iter);		
}




/** ComputeStateActionValues
   
    threshold - exit when difference in Q is smaller than the threshold
    max_iter - exit when the number of iterations reaches max_iter

*/

void MultiMDPValueIteration::ComputeStateActionValues(real threshold, int max_iter)
{
#if 1
    ComputeStateValues(threshold, max_iter);
#else
    // this implementation fails!
    action_counts.Resize(n_states, n_actions);
    action_counts.Clear();

    int n_iter = 0;
    do {
        pV_xi = V_xi;
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
                Q_xi(s,a) = Q_sa;
            }
            V_xi(s) = Q_a_max;
            //printf ("%d ", a_max[s]);
            action_counts(s, a_max[s]) += 1.0;
            //printf ("%f ", V_xi(s));
        }
        //printf ("# opt_act\n");

        // Calculate the value of each MDP for the current best average action
        for (int i=0; i<n_mdps; ++i) {
            Vector tmpV = V[i];
            for (int s=0; s<n_states; s++) {
                tmpV(s) = ComputeStateActionValueForSingleMDP(i, s, a_max[s]);
            }
            V[i] = tmpV;
        }
        if (max_iter > 0) {
            max_iter--;
        }
        Delta = abs(V_xi - pV_xi).Sum();
        n_iter++;
		//printf("%f # delta\n", Delta);
        //V.print(stdout);
        //pV.print(stdout);
    } while(Delta >= threshold && max_iter != 0);
    logmsg("Exiting at delta :%f, iter :%d\n", Delta, n_iter);	
#endif
}

/// Get the policy that is optimal for the mixed MDP.
FixedDiscretePolicy* MultiMDPValueIteration::getPolicy()
{
    FixedDiscretePolicy* policy = new FixedDiscretePolicy(n_states, n_actions);
#if 1
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
#else
    for (int s=0; s<n_states; s++) {
        Vector c = action_counts.getRow(s);
        c /= c.Sum();
        Vector* p = policy->getActionProbabilitiesPtr(s);
        for (int a=0; a<n_actions; a++) { 
            (*p)[a] = c(a);
        }
    }
#endif
    return policy;
}


