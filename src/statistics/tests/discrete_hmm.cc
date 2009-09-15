#include "DiscreteHiddenMarkovModel.h"
#include "DiscreteHiddenMarkovModelPF.h"
#include "Matrix.h"
#include "Dirichlet.h"
#include "Random.h"

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

struct Stats
{
    std::vector<real> loss;
    std::vector<real> accuracy;
    Stats(int T) : loss(T), accuracy(T)
    {
        for (int t=0; t<T; ++t) {
            loss[t] = 0;
            accuracy[t] = 0;
        }
    }
};

void TestBelief (DiscreteHiddenMarkovModel* hmm,
                 int T,
                 Stats& state_stats,
                 Stats& observation_stats,
                 Stats& pf_stats,
                 Stats& pf_mix_stats)
{
    real threshold = 0.5;
    real stationarity = 0.9;
    int n_particles = 1024;
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
            observation_stats.loss[t] += 1;
        }
        observation_stats.accuracy[t] += Px_t[x];

        // adapt belief state to observation
        hmm_belief_state.Observe(x);

        // see if current state is tracked
        Vector Ps_t = hmm_belief_state.getBelief();
        int predicted_state = ArgMax(Ps_t);

        if (predicted_state!=s) {
            state_stats.loss[t] += 1;
        }
        state_stats.accuracy[t] += Ps_t[s];
        
        // --- particle filter ---
        // predict observation by PF
        Px_t = hmm_pf.getPrediction();
        predicted_observation = ArgMax(Px_t);
        if (predicted_observation != x) {
            pf_stats.loss[t] += 1;
        }
        pf_stats.accuracy[t] += Px_t[x];

        // adapt PF
        hmm_pf.Observe(x);

        // --- particle filter mixture ---
        // predict observation by PFM
        Px_t = hmm_pf_mixture.getPrediction();
        predicted_observation = ArgMax(Px_t);
        if (predicted_observation != x) {
            pf_mix_stats.loss[t] += 1;
        }
        pf_mix_stats.accuracy[t] += Px_t[x];

        // adapt PF
        hmm_pf_mixture.Observe(x);

    }

    printf("Real model\n");
    hmm->Show();
}


int main(int argc, char** argv)
{
    if (argc != 4) {
        fprintf(stderr, "Usage: discrete_hmm n_states n_observations stationarity\n");
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
    
    Vector x(10);
    Vector y = exp(x);

    int T = 1000;
    Stats state_stats(T);
    Stats observation_stats(T);
    Stats pf_stats(T);
    Stats pf_mix_stats(T);
    
    int n_iter = 1;
    for (int i=0; i<n_iter; ++i) {
        DiscreteHiddenMarkovModel* hmm = MakeRandomDiscreteHMM(n_states,  n_observations, stationarity);
        TestBelief(hmm, T, state_stats, observation_stats, pf_stats, pf_mix_stats);
        delete hmm;
    }
    
    real norm = 1.0 / (real) n_iter;

    for (int t=0; t<T; ++t) {
        printf ("%f %f %f %f %f %f %f %f # loss\n", 
                state_stats.loss[t]*norm,
                state_stats.accuracy[t]*norm,
                observation_stats.loss[t]*norm,
                observation_stats.accuracy[t]*norm,
                pf_stats.loss[t]*norm,
                pf_stats.accuracy[t]*norm,
                pf_mix_stats.loss[t]*norm,
                pf_mix_stats.accuracy[t]*norm);

    }
    return 0;
}

