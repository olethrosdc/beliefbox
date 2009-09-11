#include "DiscreteHiddenMarkovModel.h"
#include "Matrix.h"
#include "Dirichlet.h"
#include "Random.h"

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

    DiscreteHiddenMarkovModel hmm(Pr_S, Pr_X);
    DiscreteHiddenMarkovModelStateBelief hmm_belief_state(n_states);
    hmm_belief_state.hmm = &hmm;

    int T = 100;

    for (int t=0; t<T; ++t) {
        int x = hmm.generate();
        int s = hmm.getCurrentState();
        hmm_belief_state.Observe(x);
        Vector Pt = hmm_belief_state.getBelief();
        printf ("%d %d ", x, s);
        printf("%f ", Pt[s]);
        for (int i=0; i<n_states; ++i) {
            printf("%.4f ", Pt[i]);
        }
        printf("\n");
    }
    return 0;
}
