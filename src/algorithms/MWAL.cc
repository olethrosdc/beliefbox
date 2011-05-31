// -*- Mode: C++ -*-

// Multiplicative weights for apprenticeship learning (MWAL) algorithm.

#include "MWAL.h"
#include "DiscretePolicy.h"
#include "ValueIteration.h"

/// Calculate the feature counts from a given set of demonstrations.
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

/// Calculate the feature expectations for a particular policy and
/// discount factor gamma, up to accuracy epsilon.
Vector MWAL::CalculateFeatureExpectation(DiscreteMDP& mdp,
                                         FixedDiscretePolicy& policy,
                                         real gamma, real epsilon)
{
    // Calculate Markov chain transitions
    // P(i,j) = Pr(s_{t+1} = j | s_t = i).
    Matrix P(n_states, n_states); 
    P.Clear();
    for (int i=0; i<n_states; ++i) {
        Vector Pi = policy.getActionProbabilities(i);
        for (int a=0; a<n_actions; ++a) {
            real P_a = Pi(a);
            for (int j=0; j<n_states; ++j) {
                real p = mdp.getTransitionProbability(i, a, j);
                P(j, i) += P_a * p;
            }
        }
    }
    const Matrix& Pc = P;    
    // Calculate mu If D is a Nx1 matrix, P is a N*N matrix, then PD
    // is N x 1. If each column of P is a probability distribution
    // 
    Vector mu(n_states);
    int T = (int) (log(epsilon * (1.0 - gamma)) / log(gamma));
    real discount = 1;
    Vector D(Vector::Unity(n_states));
    D /= (real) n_states;
    for (int t=0; t<T; ++t) {
        mu += D * discount;
        const Vector& Dc = D;
        D = Pc * Dc;
        discount *= gamma;
    }
    
    return mu;
}


/// Compute a policy for a particular mdp, with discount factor gamma,
/// an error tolerance and a particular number of iterations T.
void MWAL::Compute(DiscreteMDP& original_mdp, real gamma, real epsilon, int T)
{
    DiscreteMDP& mdp(original_mdp); // make a copy
    // Setup
    real beta = 1.0 / (1.0 + sqrt(2.0 * log((real) n_states) / (real) T));
    Vector W = Vector::Unity(n_states);
    std::vector<Distribution*> rewards(n_states);
    

    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            mean_policy.p[s](a) = 0;
        }
    }

    
    //printf ("MU_E: "); mu_E.print(stdout);
    // main loop
    if (T < 0) {
        T = (int) (1.0 / (epsilon * (1.0 - gamma)));
        T *= T;
    }
    for (int t=0; t<T; ++t) {
        Vector w = W / W.Sum();
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                mdp.setFixedReward(s, a, w(s));
            }
        }
        
        ValueIteration VI (&mdp, gamma);
        VI.ComputeStateActionValues(epsilon, 100);
        
        FixedDiscretePolicy policy(n_states, n_actions, VI.Q);
        Vector mu_t = CalculateFeatureExpectation(mdp, policy, gamma, epsilon);
        Vector G = ((mu_t - mu_E) * (1.0 - gamma) + 2.0) / 4;
        W *= exp(G * log(beta));

        real delta = 0;
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                real p_sa = policy.getActionProbability(s, a);
                delta += fabs(mean_policy.p[s](a) / (real) (1 + t) - p_sa);
                mean_policy.p[s](a) += p_sa;
            }
        }
        if (delta < epsilon ) {
            printf ("Breaking at %d instead of %d\n", t, T);
            T = t;
        }
        //printf ("delta: %f\n", delta);
        //printf ("MU_t: "); mu_t.print(stdout);
        //policy.Show();
    }
    real invT = 1.0 / (real) T;
    
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a) {
            mean_policy.p[s](a) *= invT;
        }
    }
    
    //mean_policy.Show();
    for (int s=0; s<n_states; ++s) {
        delete rewards[s];
    }
}
