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

#include "DiscreteHiddenMarkovModelPF.h"
#include "Dirichlet.h"
#include "Matrix.h"
#include <cmath>

#undef REPLACE_ALL_LOW

DiscreteHiddenMarkovModelPF::DiscreteHiddenMarkovModelPF(real threshold, real stationarity, int n_states_, int n_observations_, int n_particles_)
    :  n_states(n_states_), n_observations(n_observations_), 
       n_particles(n_particles_), hmm(n_particles), belief(n_particles),
       P_x(n_particles), log_P_x(n_particles),
       w(n_particles), log_w(n_particles),
       state_prior(n_states),
       observation_prior(n_states)
{

    //printf ("%d states\n", n_states);
    //    state_prior.resize(n_states);
    //    observation_prior.resize(n_states);
    
    // initialise the prior matrices
    for (int i=0; i<n_states; ++i) {
        state_prior[i] = new DirichletDistribution(n_states, threshold);
        state_prior[i]->Alpha(i) = threshold * stationarity * (real) (n_states - 1) / (1.0 - stationarity);
        observation_prior[i] = new DirichletDistribution(n_observations, threshold);
    }
    
    // initialise the particle weights
    real prior = 1.0 / (real) n_particles;
    real log_prior = log (prior);
    for (int k=0; k<n_particles; ++k) {
        w[k] = prior;
        log_w[k] = log_prior;
    }

    // initialise the particles
    for (int k=0; k<n_particles; ++k) {
        // initialise state transition matrix
        Matrix P_S(n_states, n_states);
        for (int i=0; i<n_states; ++i) {
            Vector p = state_prior[i]->generate();
            for (int j=0; j<n_states; ++j) {
                P_S(i,j) = p[j];
            }
        }
        // initialise observation matrix
        Matrix P_X(n_states, n_observations);
        for (int i=0; i<n_states; ++i) {
            Vector p = observation_prior[i]->generate();
            for (int j=0; j<n_observations; ++j) {
                P_X(i,j) = p[j];
            }
        }
        // create HMM
        hmm[k] = new DiscreteHiddenMarkovModel(P_S, P_X);
        belief[k] = new DiscreteHiddenMarkovModelStateBelief(hmm[k]);
    }
}

DiscreteHiddenMarkovModelPF::~DiscreteHiddenMarkovModelPF()
{
    //int max_k = ArgMax(w);
    //hmm[max_k]->Show();
    for (int i=0; i<n_states; ++i) {
        delete state_prior[i];
        delete observation_prior[i];
    }
    for (int k=0; k<n_particles; ++k) {
        delete hmm[k];
        delete belief[k];
    }
}


real DiscreteHiddenMarkovModelPF::Observe(int x)
{
    real log_sum = LOG_ZERO;
    
    // calculate p(x|k) and p(x) = sum_k p(x,k)
   
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);

    return exp(log_sum);
}

Vector DiscreteHiddenMarkovModelPF::getPrediction()
{
    Vector p_x(n_observations);
    for (int i=0; i<n_observations; ++i) {
        p_x[i] = 0.0;
    }
    // p(k|x) = p(x|k) / p(x)
    for (int k=0; k<n_particles; ++k) {
        p_x += belief[k]->getPrediction() * w[k];
    }

    return p_x;
}

void DiscreteHiddenMarkovModelPF::Reset()
{
}

void DiscreteHiddenMarkovModelPF::Show()
{
    for (int k=0; k<n_particles; ++k) {
        printf ("w[%d] = %f\n", k, w[k]);
        hmm[k]->Show();
    }
}

//------------------ DiscreteHiddenMarkovModelPF_ReplaceLowest ---------------//


real DiscreteHiddenMarkovModelPF_ReplaceLowest::Observe(int x)
{
    int min_k = ArgMin(log_w);
#ifdef REPLACE_ALL_LOW
    int reps = n_particles;
    while (log_w[min_k] < replacement_threshold && reps-- > 0) {
#else
    if (log_w[min_k] < replacement_threshold) {
#endif
        int k = DiscreteDistribution::generate(w);
        real alpha = 0.1;
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        std::vector<MultinomialDistribution>& PS_min = hmm[min_k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_min = hmm[min_k]->getObservationProbablities();
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_states; ++j) {
                PS_min[i].Pr(j) = PS_k[i].Pr(j) * alpha + PS_min[i].Pr(j) * (1 - alpha);
            }
            for (int j=0; j<n_observations; ++j) {
                PX_min[i].Pr(j) = PX_k[i].Pr(j) * alpha + PX_min[i].Pr(j) * (1 - alpha);
            }
        }
        log_w[min_k] = logAdd(log(alpha) + log_w[k], log(1 - alpha) + log_w[min_k]);
        min_k = ArgMin(log_w);
    }
    // normalise weights
    log_w -= log_w.logSum();
        
    real log_sum = LOG_ZERO;

    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);


    return exp(log_sum);
}

//------------------ DiscreteHiddenMarkovModelPF_ISReplaceLowest ---------------//

/** Replace particles under threshold, but adjust the weight first.
    
    This is the same as the standard filter.

    First, calculate
    \f[
    P_{t+1}(q_i) = P_t (q_i | x_{t+1}).
    \f]
    Then replace low-weight particles by doing
 */
real DiscreteHiddenMarkovModelPF_ISReplaceLowest::Observe(int x)
{
    real log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;

    int min_k = ArgMin(log_w);
#ifdef REPLACE_ALL_LOW
    int reps = n_particles;
    while (log_w[min_k] < replacement_threshold && reps-- > 0) {
#else
    if (log_w[min_k] < replacement_threshold) {
#endif
            /// mix weight with k
        int k = DiscreteDistribution::generate(w);
        real alpha = 0.1;
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        std::vector<MultinomialDistribution>& PS_min = hmm[min_k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_min = hmm[min_k]->getObservationProbablities();
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_states; ++j) {
                PS_min[i].Pr(j) = PS_k[i].Pr(j) * alpha + PS_min[i].Pr(j) * (1 - alpha);
            }
            for (int j=0; j<n_observations; ++j) {
                PX_min[i].Pr(j) = PX_k[i].Pr(j) * alpha + PX_min[i].Pr(j) * (1 - alpha);
            }
        }
        log_w[min_k] = logAdd(log(alpha) + log_w[k], log(1 - alpha) + log_w[min_k]);
        min_k = ArgMin(log_w);
    }
    // normalise weights
    log_w -= log_w.logSum();
        
    log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);


    return exp(log_sum);
}



//------------ DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichlet ---------//

/** Replace particles under threshold, but adjust the weight first.
    
    This is the same as the standard filter.

    First, calculate
    \f[
    P_{t+1}(q_i) = P_t (q_i | x_{t+1}).
    \f]
    Then replace low-weight particles by sampling from the Dirichlet
    mixture.
 */
real DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichlet::Observe(int x)
{
    T++;
    real scale = (real) T;
    real log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;

    int min_k = ArgMin(log_w);
#ifdef REPLACE_ALL_LOW
    int reps = n_particles;
    while (log_w[min_k] < replacement_threshold && reps-- > 0) {
#else
    if (log_w[min_k] < replacement_threshold) {
#endif
            /// mix weight with k
        int k = DiscreteDistribution::generate(w);
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        std::vector<MultinomialDistribution>& PS_min = hmm[min_k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_min = hmm[min_k]->getObservationProbablities();
            // create Dirichlet and sample
        for (int i=0; i<n_states; ++i) {
                // state Dirichlet
            Vector v_s(n_states);
            for (int j=0; j<n_states; ++j) {
                v_s[j] = scale * PS_k[i].Pr(j);
            }
            DirichletDistribution dir_s(v_s);
            dir_s.generate(v_s);
            for (int j=0; j<n_states; ++j) {
                PS_min[i].Pr(j) = v_s[j];
            }
            
                // observation Dirichlet
            Vector v_x(n_observations);
            for (int j=0; j<n_observations; ++j) {
                v_x[j] = scale * PX_k[i].Pr(j);
            }
            DirichletDistribution dir_x(v_x);
            dir_x.generate(v_x);
            for (int j=0; j<n_observations; ++j) {
                PX_min[i].Pr(j) = v_x[j];
            }
        }
        log_w[min_k] = log_w[k] - (real) n_particles;
        min_k = ArgMin(log_w);
    }
    // normalise weights
    log_w -= log_w.logSum();
        
    log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);


    return exp(log_sum);
}


//----- DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichletExact -----------//
real DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichletExact::Observe(int x)
{
    T++;
    history.push_back(x);
    real scale = (real) T;
    real log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;

    int min_k = ArgMin(log_w);
#ifdef REPLACE_ALL_LOW
    int reps = n_particles;
    while (log_w[min_k] < replacement_threshold && reps-- > 0) {
#else
    if (log_w[min_k] < replacement_threshold) {
#endif
            /// mix weight with k
        int k = DiscreteDistribution::generate(w);
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        std::vector<MultinomialDistribution>& PS_min = hmm[min_k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_min = hmm[min_k]->getObservationProbablities();
            // create Dirichlet and sample
        for (int i=0; i<n_states; ++i) {
                // state Dirichlet
            Vector v_s(n_states);
            for (int j=0; j<n_states; ++j) {
                v_s[j] = scale * PS_k[i].Pr(j);
            }
            DirichletDistribution dir_s(v_s);
            dir_s.generate(v_s);
            for (int j=0; j<n_states; ++j) {
                PS_min[i].Pr(j) = v_s[j];
            }
            
                // observation Dirichlet
            Vector v_x(n_observations);
            for (int j=0; j<n_observations; ++j) {
                v_x[j] = scale * PX_k[i].Pr(j);
            }
            DirichletDistribution dir_x(v_x);
            dir_x.generate(v_x);
            for (int j=0; j<n_observations; ++j) {
                PX_min[i].Pr(j) = v_x[j];
            }
        }
        log_w[min_k] = log_w[k] - (real) n_particles;
        belief[min_k]->Reset();
        belief[min_k]->Observe(history);
        min_k = ArgMin(log_w);
    }
    // normalise weights
    log_w -= log_w.logSum();
        
    log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);


    return exp(log_sum);
}


//----- DiscreteHiddenMarkovModelPF_ISReplaceLowestDirichletExact -----------//
real DiscreteHiddenMarkovModelPF_BootstrapDirichletExact::Observe(int x)
{
    // TODO
}


//------------- DiscreteHiddenMarkovModelPF_ReplaceLowestExact --------------//
real DiscreteHiddenMarkovModelPF_ReplaceLowestExact::Observe(int x)
{
    history.push_back(x);
    int min_k = ArgMin(log_w);
#ifdef REPLACE_ALL_LOW
    int reps = n_particles;
    while (log_w[min_k] < replacement_threshold && reps-- > 0) {
#else
    if (log_w[min_k] < replacement_threshold) {
#endif
        int k = DiscreteDistribution::generate(w);
        real alpha = 0.1;
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        std::vector<MultinomialDistribution>& PS_min = hmm[min_k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_min = hmm[min_k]->getObservationProbablities();
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_states; ++j) {
                PS_min[i].Pr(j) = PS_k[i].Pr(j) * alpha + PS_min[i].Pr(j) * (1 - alpha);
            }
            for (int j=0; j<n_observations; ++j) {
                PX_min[i].Pr(j) = PX_k[i].Pr(j) * alpha + PX_min[i].Pr(j) * (1 - alpha);
            }
        }
        log_w[min_k] = logAdd(log(alpha) + log_w[k], log(1 - alpha) + log_w[min_k]);
        belief[min_k]->Reset();
        belief[min_k]->Observe(history);
        min_k = ArgMin(log_w);
    }
    // normalise weights
    log_w -= log_w.logSum();
   
    real log_sum = LOG_ZERO;

    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);

    return exp(log_sum);
}

//------------ DiscreteHiddenMarkovModelPF_ISReplaceLowestExact -------------//
// Importance Sampling
// Lowest replacement
// Exact state belief update
real DiscreteHiddenMarkovModelPF_ISReplaceLowestExact::Observe(int x)
{

   real log_sum = LOG_ZERO;
    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;


    history.push_back(x);
    int min_k = ArgMin(log_w);
#ifdef REPLACE_ALL_LOW
    int reps = n_particles;
    while (log_w[min_k] < replacement_threshold && reps-- > 0) {
#else
    if (log_w[min_k] < replacement_threshold) {
#endif
        int k = DiscreteDistribution::generate(w);
        real alpha = 0.1;
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        std::vector<MultinomialDistribution>& PS_min = hmm[min_k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_min = hmm[min_k]->getObservationProbablities();
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_states; ++j) {
                PS_min[i].Pr(j) = PS_k[i].Pr(j) * alpha + PS_min[i].Pr(j) * (1 - alpha);
            }
            for (int j=0; j<n_observations; ++j) {
                PX_min[i].Pr(j) = PX_k[i].Pr(j) * alpha + PX_min[i].Pr(j) * (1 - alpha);
            }
        }
        log_w[min_k] = logAdd(log(alpha) + log_w[k], log(1 - alpha) + log_w[min_k]);
        belief[min_k]->Reset();
        belief[min_k]->Observe(history);
        min_k = ArgMin(log_w);
    }
    // normalise weights
    log_w -= log_w.logSum();
   
    log_sum = LOG_ZERO;

    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);

    return exp(log_sum);
}




//------------------ DiscreteHiddenMarkovModelPF_Resample ---------------//



/// Here we resample from the Dirichlet created by the PF.
real DiscreteHiddenMarkovModelRBPF::Observe(int x)
{
    t++; // increase number of observations

    // set up distribution to sample from
    std::vector<DirichletDistribution> dS(n_states);
    std::vector<DirichletDistribution> dX(n_states);

    // fill in sampling distribution values
    for (int k=0; k<n_particles; ++k) {
        //real alpha = sqrt(real (t));
        std::vector<MultinomialDistribution>& PS_k = hmm[k]->getStateProbablities();
        std::vector<MultinomialDistribution>& PX_k = hmm[k]->getObservationProbablities();
        for (int i=0; i<n_states; ++i) {
            for (int j=0; j<n_states; ++j) {
                dS[i].Alpha(j) += PS_k[i].Pr(j);
            }
            for (int j=0; j<n_observations; ++j) {
                dX[i].Alpha(j) += PX_k[i].Pr(j);
            }
        }
    }
        
    real log_sum = LOG_ZERO;

    // calculate p(x|k) and p(x) = sum_k p(x,k)
    for (int k=0; k<n_particles; ++k) {
        P_x[k] = belief[k]->Observe(x);
        log_P_x[k] = log(P_x[k]) + log_w[k];
        log_sum = logAdd(log_sum, log_P_x[k]);
    }
    
    // p(k|x) = p(x|k) / p(x)
    log_w = log_P_x - log_sum;
    w = exp(log_w);


    return exp(log_sum);
}

void DiscreteHiddenMarkovModelRBPF::Reset()
{
    t = 0;
}
