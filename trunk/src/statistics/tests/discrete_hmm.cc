#include "DiscreteHiddenMarkovModel.h"
#include "DiscreteHiddenMarkovModelPF.h"
#include "Matrix.h"
#include "Dirichlet.h"
#include "Random.h"
#include "CumulativeStats.h"

DiscreteHiddenMarkovModel* MakeRandomDiscreteHMM(int n_states, int n_observations, real stationarity)
{
    assert (n_states > 0);
    assert (n_observations > 0);
    assert (stationarity >= 0 && stationarity <= 1);

    Matrix Pr_S(n_states, n_states);
    Matrix Pr_X(n_states, n_observations);
    for (int i=0; i<n_states; ++i) {
        real sum = 0.0;
        for (int j=0; j<n_observations; ++j) {
            Pr_X(i,j) = 0.1*true_random(false);
            if (i==j) {
                Pr_X(i,j) += 1.0;
            }
            sum += Pr_X(i,j);
        }
        for (int j=0; j<n_observations; ++j) {
            Pr_X(i,j) /=  sum;
        }
        
    }
    Matrix S(n_states, n_states);
    for (int src=0; src<n_states; ++src) {
        Vector P(n_states);
        for (int i=0; i<n_states; ++i) {
            if (i<=src) {
                P[i] = exp(true_random(false));
            } else {
                P[i] = exp(10.0*true_random(false));
            }
        }
        P /= P.Sum();
        P *= (1 - stationarity);
        P[src] += stationarity;
        P /= P.Sum();
            //real sum = 0.0;
        for (int dst=0; dst<n_states; ++dst) {
            Pr_S(src,dst) = P[dst];
        }
    }

    return new DiscreteHiddenMarkovModel (Pr_S, Pr_X);
}



void TestBelief (DiscreteHiddenMarkovModel* hmm,
                 int T,
                 real threshold,
                 real stationarity,
                 int n_particles,
                 CumulativeStats& state_stats,
                 CumulativeStats& observation_stats,
                 CumulativeStats& pf_stats,
                 CumulativeStats& pf_mix_stats)
{
    DiscreteHiddenMarkovModelStateBelief hmm_belief_state(hmm);
    DiscreteHiddenMarkovModelPF hmm_pf(threshold, stationarity, hmm->getNStates(), hmm->getNObservations(), n_particles);

    int max_states = 8; //Max(16, 2 * hmm->getNStates())

    DHMM_PF_Mixture hmm_pf_mixture(threshold, stationarity, hmm->getNObservations(), n_particles, max_states);

    for (int t=0; t<T; ++t) {
        // perdict next observation
        Vector Px_t = hmm_belief_state.getPrediction();
        int predicted_observation = ArgMax(Px_t);

        // generate next observation and get state
        int x = hmm->generate();
        int s = hmm->getCurrentState();
        
        // add observation error
        if (predicted_observation != x) {
            observation_stats.SetValue(t, 1);
        } else {
            observation_stats.SetValue(t, 0);
        }

        // adapt belief state to observation
        hmm_belief_state.Observe(x);

        // see if current state is tracked
        Vector Ps_t = hmm_belief_state.getBelief();
        int predicted_state = ArgMax(Ps_t);

        if (predicted_state!=s) {
            state_stats.SetValue(t, 1);
        } else {
            state_stats.SetValue(t, 0);
        }

        
        // --- particle filter ---
        // predict observation by PF
        Px_t = hmm_pf.getPrediction();
        predicted_observation = ArgMax(Px_t);
        if (predicted_observation != x) {
            pf_stats.SetValue(t, 1);
        } else {
            pf_stats.SetValue(t, 0);
        }

        // adapt PF
        hmm_pf.Observe(x);

        // --- particle filter mixture ---
        // predict observation by PFM
        Px_t = hmm_pf_mixture.getPrediction();
        predicted_observation = ArgMax(Px_t);
        if (predicted_observation != x) {
            pf_mix_stats.SetValue(t, 1);
        } else {
            pf_mix_stats.SetValue(t, 0);
        }

        // adapt PF
        hmm_pf_mixture.Observe(x);

    }

    //printf("Real model\n");
    //hmm->Show();
}


int main(int argc, char** argv)
{
    if (argc != 7) {
        fprintf(stderr, "Usage: discrete_hmm n_states n_observations stationarity n_particles T n_iter\n");
        return -1;
    }
    int n_states = atoi(argv[1]);
    if (n_states <= 0) {
        fprintf (stderr, "Invalid number of states %d\n", n_states);
    }

    int n_observations = atoi(argv[2]);
    if (n_observations <= 0) {
        fprintf (stderr, "Invalid number of states %d\n", n_observations);
    }

    real stationarity = atof(argv[3]);
    if (stationarity < 0 || stationarity > 1) {
        fprintf (stderr, "Invalid stationarity %f\n", stationarity);
    }

    int n_particles = atoi(argv[4]);
    if (n_particles <= 0) {
        fprintf (stderr, "Invalid n_particles %d\n", n_particles);
    }

    int T = atoi(argv[5]);
    if (T <= 0) {
        fprintf (stderr, "Invalid T %d\n", T);
    }

    int n_iter = atoi(argv[6]);
    if (n_iter <= 0) {
        fprintf (stderr, "Invalid n_iter %d\n", n_iter);
    }
    

    real threshold = 0.5; // threshold for the prior in the estimated HMMs.

    Vector x(10);
    Vector y = exp(x);


    CumulativeStats state_stats(T, n_iter);
    CumulativeStats observation_stats(T, n_iter);
    CumulativeStats pf_stats(T, n_iter);
    CumulativeStats pf_mix_stats(T, n_iter);
    
    for (int i=0; i<n_iter; ++i) {
        state_stats.SetSequence(i);
        observation_stats.SetSequence(i);
        pf_stats.SetSequence(i);
        pf_mix_stats.SetSequence(i);
        fprintf (stderr, "Iter: %d / %d\n", i + 1, n_iter);
        DiscreteHiddenMarkovModel* hmm = MakeRandomDiscreteHMM(n_states,  n_observations, stationarity);
        TestBelief(hmm, T, threshold, stationarity, n_particles, state_stats, observation_stats, pf_stats, pf_mix_stats);
        delete hmm;
    }
    
    
    real percentile = 0.1;
    
    //Matrix M = pf_stats.C;
    CumulativeStats pf_to_mix_stats(pf_stats.C - pf_mix_stats.C);

    observation_stats.Sort();
    Vector hmm_mean = observation_stats.Mean();
    Vector hmm_top = observation_stats.TopPercentile(percentile);
    Vector hmm_bottom = observation_stats.BottomPercentile(percentile); 
    
    pf_stats.Sort();
    Vector pf_mean = pf_stats.Mean();
    Vector pf_top = pf_stats.TopPercentile(percentile);
    Vector pf_bottom = pf_stats.BottomPercentile(percentile);
    
    pf_mix_stats.Sort();
    Vector pf_mix_mean = pf_mix_stats.Mean();
    Vector pf_mix_top = pf_mix_stats.TopPercentile(percentile);
    Vector pf_mix_bottom = pf_mix_stats.BottomPercentile(percentile);
    for (int t=0; t<T; ++t) {
        printf ("%f %f %f %f %f %f %f %f %f # loss\n", 
                hmm_bottom[t], hmm_mean[t], hmm_top[t],
                pf_bottom[t], pf_mean[t], pf_top[t],
                pf_mix_bottom[t], pf_mix_mean[t], pf_mix_top[t]);
    }
    

    pf_to_mix_stats.Sort();
    Vector pf_to_mix_mean = pf_to_mix_stats.Mean();
    Vector pf_to_mix_top = pf_to_mix_stats.TopPercentile(percentile);
    Vector pf_to_mix_bottom = pf_to_mix_stats.BottomPercentile(percentile);
    for (int t=0; t<T; ++t) {
        printf ("%f %f %f # pf_to_mix\n", 
                pf_to_mix_bottom[t], pf_to_mix_mean[t], pf_to_mix_top[t]);
    }

    return 0;
}

