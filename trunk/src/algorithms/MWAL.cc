// -*- Mode: C++ -*-

// Multiplicative weights for apprenticeship learning (MWAL) algorithm.

#include "MWAL.h"
#include "DiscretePolicy.h"
#include "ValueIteration.h"

void MWAL::CalculateFeatureCounts(Demonstrations<int, int>& D)
{

    for (int s=0; s<n_states; ++s) {
        mu_E(s) = 0;
    }


    int K = D.trajectories.size();
    for (int k=0; k<K; ++k) {
        real discount = 1;
        for (uint t=0; t<D.trajectories[k].x.size(); ++t) {
            mu_E(D.trajectories[k].x[t].first) += discount;
            discount *= gamma;
        }
    }
    mu_E /= (real) K;
}

Vector MWAL::CalculateFeatureExpectation(DiscreteMDP& mdp,
                                         FixedDiscretePolicy& policy,
                                         real gamma, real epsilon)
{
    // Calculate Markov chain transitions
    // P(j,i) = Pr(s_{t+1} = j | s_t = i).
    Matrix P(n_states, n_states); 
    for (int i=0; i<n_states; ++i) {
        Vector Pi = policy.getActionProbabilities(i);
        for (int a=0; a<n_actions; ++a) {
            real P_a = Pi(a);
            for (int j=0; j<n_states; ++j) {
                P(j, i) += P_a * mdp.getTransitionProbability(i, a, j);
            }
        }
    }

    const Matrix& Pc = P;    
    // Calculate mu
    Vector mu(n_states);
    int T = (int) (log(epsilon) / log(gamma));
    real discount = 1;
    Vector D(Vector::Unity(n_states));
    D /= (real) n_states;
    for (int t=0; t<T; ++T) {
        mu += D * discount;
        const Vector& Dc = D;
        D = Pc * Dc;
        discount *= gamma;
    }
    
    return mu;
}


void MWAL::Compute(DiscreteMDP& mdp, real gamma, real epsilon, real T)
{
    // Setup
    real beta = 1.0 / (1.0 + sqrt(2.0 * log((real) n_states) / (real) T));
    Vector W(Vector::Unity(n_states));
    std::vector<Distribution*> rewards(n_states);
    for (int s=0; s<n_states; ++s) {
        rewards[s] = new SingularDistribution(0);
        for (int a=0; a<n_actions; ++a) {
            mdp.setRewardDistribution(s, a, rewards[s]);
        }
    }
    
    
    FixedDiscretePolicy mean_policy(n_states, n_actions);
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            mean_policy.p[s](a) = 0;
        }
    }

    // main loop
    for (int t=0; t<T; ++t) {
        Vector w = W / W.Sum();
        for (int s=0; s<n_states; ++s) {
            rewards[s]->setMean(w(s));
        }
        
        ValueIteration VI (&mdp, gamma, epsilon);
        FixedDiscretePolicy policy(n_states, n_actions, VI.Q);
        Vector mu_t = CalculateFeatureExpectation(mdp, policy, gamma, epsilon);
        Vector G = ((mu_t - mu_E) * (1.0 - gamma) + 2.0) / 4;
        W *= exp(G * log(beta));
        
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                mean_policy.p[s](a) += policy.getActionProbability(s, a);
            }
        }
    }
    real invT = 1.0 / (real) T;
    
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            mean_policy.p[s](a) *= invT;
        }
    }
    
    for (int s=0; s<n_states; ++s) {
        delete rewards[s];
    }


    
}
