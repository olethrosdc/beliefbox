// -*- Mode: c++ -*-
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: test_value_iteration.c,v 1.8 2006/11/12 23:16:03 olethros Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "DiscreteMDP.h"
#include "ValueIteration.h"
#include "Distribution.h"
#include "Vector.h"
#include "EasyClock.h"
#include "SingularDistribution.h"
#include <ctime>
#include <cstdlib>

class SimpleBelief
{
public:
    real y[2];
    real w[2];
    real w2[2];
    int N;
    SimpleBelief(real ya, real yb, real wa, real wb)
    {
        y[0] = ya;
        y[1] = yb;
        w[0] = wa;
        w[1] = wb;
        this->N = 2;
        //reward_distribution.p = getMean();
    }
    void Observe(real x)
    {
        // Evaluate the new posterior up to a normalising constant
        real sum = 0.0f;
        for (int i=0; i<N; i++) {
            real q = y[i];
            real likelihood = q*x + (1-q)*(1-x);
            real prior = w[i];
            w2[i] = likelihood*prior;
            sum += w2[i];
        }
        if (sum==0.0f) {
            fprintf (stderr, "ER! ERROR: 0 mass on prior %d!!\n", N);
            for (int i=0; i<N; i++) {
                fprintf (stderr, "%d -> %f, %f -> %f\n", i, y[i], w[i], w2[i]);
            }
            fflush(stderr);
            fprintf (stderr, "ER! ERROR: 0 mass on prior!!\n");
            exit(-1);
        }

        // Normalise to create a new filtering distribution
        real isum = 1.0f / sum;
        for (int i=0; i<N; i++) {
            w[i] = w2[i] * isum;
        }
        //		reward_distribution.p = getMean();
    }
    real getMean()
    {
        real mean = 0.0f;
        for (int i=0; i<N; i++) {
            mean += w[i] * y[i];
        }
        return mean;
    }
};

int test_grid_world();
int test_belief_mdp(int T, real gamma, real prior);
int play_belief_mdp(int T, real gamma, SimpleBelief& initial_belief);
real evaluate_belief_mdp(int T, real gamma, int n_trials, int horizon);
real test_greedy(int n_trials, int horizon, real alpha);
real test_hoeffding_bound(int n_trials, int horizon, real delta);

int main(int argc, char** argv)
{
    srand48(time(NULL));
    srand(time(NULL));
    if (argc<2) {
      printf ("args:\n"
              "test_belief T gamma prior\n"
	      "belief T gamma n_trials horizon\n"
	      "grid\n"
              "greedy n_trials horizon alpha\n"
	      "hoeffding n_trials horizon delta\n");
      exit(-1);
    }

    if (!strcmp("test_belief", argv[1])) {
        if (argc!=5) {
            fprintf (stderr, "argc\n");
            exit(-1);
        }
        return test_belief_mdp(atoi(argv[2]), atof(argv[3]), atof(argv[4]));
    } else if (!strcmp("belief", argv[1])) {
        if (argc!=6) {
            fprintf (stderr, "argc\n");
            exit(-1);
        }
        int T = atoi(argv[2]);
        real gamma = atof(argv[3]);
        int n_trials = atoi(argv[4]);
        int horizon = atoi(argv[5]);
        printf ("# total: %f\n", evaluate_belief_mdp(T, gamma, n_trials, horizon));
        return 0;
    } else if (!strcmp("grid", argv[1])) {
        return test_grid_world();
    } else if (!strcmp("greedy", argv[1])) {
        int n_trials = atoi(argv[2]);
        int horizon = atoi(argv[3]);
        real alpha = atof(argv[4]);
        printf ("# total: %f\n", test_greedy(n_trials, horizon, alpha));
        return 0;
    } else if (!strcmp("hoeffding", argv[1])) {
        int n_trials = atoi(argv[2]);
        int horizon = atoi(argv[3]);
        real delta = atof(argv[4]);
        printf ("# total: %f\n", test_hoeffding_bound(n_trials, horizon, delta));
        return 0;
    } else {
        printf ("test_belief, belief, grid, greedy, hoeffding?\n");
        return -1;
    }
}


int test_grid_world()
{
    int n_states = 15;
    int n_actions = 4;

    DiscreteMDP mdp(n_states, n_actions, NULL, NULL);

    SingularDistribution minus_one_dist(-1.0);
    SingularDistribution zero_dist(0.5);
	
    enum Move {up=0, down, left, right};

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            mdp.setRewardDistribution(s, a, &minus_one_dist);
            for (int s2=0; s2<n_states; s2++) {
                mdp.setTransitionProbability(s, a, s2, 0.0);
            }
        }
    }

    // Set absorbing state
    mdp.setTransitionProbability(0,    up, 0, 1.0);
    mdp.setTransitionProbability(0,  down, 0, 1.0);
    mdp.setTransitionProbability(0,  left, 0, 1.0);
    mdp.setTransitionProbability(0, right, 0, 1.0);
    for (int a=0; a<n_actions; a++) {
        mdp.setRewardDistribution(0, a, &zero_dist);
    }

    // Set all other states
    for (int s=1; s<n_states; s++) {
        if (s<=4) {
            mdp.setTransitionProbability(s, up, s, 1.0);
        } else {
            mdp.setTransitionProbability(s, up, s-4, 1.0);
        }

        if (s<12) {
            if (s==11) {
                mdp.setTransitionProbability(s, down, 0, 1.0); 
            } else {
                mdp.setTransitionProbability(s, down, s+4, 1.0);
            }
        } else {
            mdp.setTransitionProbability(s, down, s, 1.0);
        }

        if (s%4) {
            mdp.setTransitionProbability(s, left, s-1, 1.0);
        } else {
            mdp.setTransitionProbability(s, left, s, 1.0);
        }
        if ((s+1)%4) {
            if (s==14) {
                mdp.setTransitionProbability(s, right, 0, 1.0);
            } else {
                mdp.setTransitionProbability(s, right, s+1, 1.0);
            }
        } else {
            mdp.setTransitionProbability(s, right, s, 1.0);
        }
    }

    real gamma = 1.0;
    ValueIteration value_iteration(&mdp, gamma, 0.5);

    value_iteration.ComputeStateValues(0.01, 100);

    printf ("---------\n");
    for (int s=0; s<n_states; s++) {
        printf ("%f", value_iteration.getValue(s));
        if ((s+1)%4) {
            printf (" ");
        } else {
            printf("\n");
        }
        fflush(stdout);
		
    }
    printf ("\n---------\n");

    value_iteration.ComputeStateActionValues(0.01, 100);

    printf ("---------\n");
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<4; a++) {
            printf ("%f ", value_iteration.getValue(s,a));
        }
        if ((s+1)%4) {
            printf ("| ");
        } else {
            printf("\n");
        }
        fflush(stdout);
		
    }
    printf ("\n---------\n");

    return 0;
}


// t = a+b, and t*(t+1)/2 always gives us the index of the
// start of the line.  Adding a gives us the complete index.
int get_bionomial_index(int a, int b)
{
    int t = a+b;
    return t*(t+1)/2 + a;
}

struct Pair
{
    int f;
    int s;
};

int test_belief_mdp(int T, real gamma, real prior)
{
    printf ("# testing belief MDP, T:%d, g:%f\n", T, gamma);
    int n_states = T*(T+1)/2;
    int n_actions = 2;

    enum Moves {STOP = 0, PLAY};

    bool add_absorbing_state = false;

    real high_value = 0.75;
    real low_value = 0.25;
    real stop_value = 0.5;
    SingularDistribution stop_reward_dist(stop_value);
    SingularDistribution play_one_dist(1.0);
    SingularDistribution play_zero_dist(0.0);

    int mdp_states = n_states;
    if (add_absorbing_state) {
        mdp_states = n_states + 1;
    }

    DiscreteMDP mdp(mdp_states, n_actions, NULL, NULL);

    ValueIteration value_iteration(&mdp, gamma, 0.5);

    std::vector<SimpleBelief*> belief(mdp_states);
    std::vector<Pair> pairs(mdp_states);
    for (int f=0; f<=T; f++) {
        for (int s=0; s<=T; s++) {
            int  t = f+s;
            if (t<T) {
                int n = get_bionomial_index(f, s);
                SMART_ASSERT(n>=0 && n<n_states)(n)(n_states);
                pairs[n].f = f;
                pairs[n].s = s;
                belief[n] = new SimpleBelief(low_value, high_value, 1.0 - prior,  prior);
                for (int i=0; i<f; i++) {
                    belief[n]->Observe(0.0);
                }
                for (int i=0; i<s; i++) {
                    belief[n]->Observe(1.0);
                }
                printf ("p(h|s=%d) = %f (prior=%f, f=%d, s=%d)\n", n, belief[n]->w[1], prior, f, s);
            }
        }
    }

    // First zero all transitions
    for (int s=0; s<mdp_states; s++) {
        for (int s2=0; s2<mdp_states; s2++) {
            mdp.setTransitionProbability(s, PLAY, s2, 0.0);
            mdp.setTransitionProbability(s, STOP, s2, 0.0);
        }
    }
    //printf ("FIRST CHECK\n");
    //mdp.Check();
	
    // When stopping, you don't get to observe nothin'
    for (int s=0; s<n_states; s++) {
        printf ("E[r|s=%d, S] = %f\n", s, stop_reward_dist.getMean());
        mdp.setRewardDistribution(s, STOP, &stop_reward_dist);
        printf ("P[s'=%d|s=%d, S] = %f\n", s, s, 1.0);
        mdp.setTransitionProbability(s, STOP, s, 1.0);
    }

    if (add_absorbing_state) {
        int s = n_states;
        printf ("E[r|s=%d, S] = %f\n", s, play_zero_dist.getMean());
        mdp.setRewardDistribution(s, STOP, &stop_reward_dist);
        printf ("P[s'=%d|s=%d, S] = %f\n", s, s, 1.0);
        mdp.setTransitionProbability(s, STOP, s, 1.0);
        printf ("E[r|s=%d, P] = %f\n", s, play_zero_dist.getMean());
        mdp.setRewardDistribution(s, PLAY, &stop_reward_dist);
        printf ("P[s'=%d|s=%d, P] = %f\n", s, s, 1.0);
        mdp.setTransitionProbability(s, PLAY, s, 1.0);
    }
    //printf ("SECOND CHECK\n");
    //mdp.Check();
    // When playing, beliefs change, until the end of time, where they are
    // stationary.
    std::vector<Distribution*> tmp_dist;
    for (int s=0; s<n_states; s++) {
        float Er = belief[s]->getMean();
        float p = belief[s]->w[1];
        Pair b=pairs[s];
        int s_success = get_bionomial_index(b.f, b.s + 1);
        int s_failure = get_bionomial_index(b.f + 1, b.s);
        BernoulliDistribution* dist = new BernoulliDistribution(Er);
        tmp_dist.push_back(dist);
        mdp.setRewardDistribution(s, PLAY, dist);
        printf ("E[r|s=%d,P] = %f\n", s, dist->getMean());
        if (s_success >= n_states || s_failure >= n_states) {
            int final_state = s;
            if (add_absorbing_state) {
                printf ("->");
                final_state = n_states;
            }
            printf ("P[s'=%d|s=%d, P] = %f\n", s, final_state, 1.0);
            mdp.setTransitionProbability(s, PLAY, final_state, 1.0);
        } else {
            printf ("P[s'=%d|s=%d, P] = %f\n", s, s_success, p);
            mdp.setTransitionProbability(s, PLAY, s_success, p);
            printf ("P[s'=%d|s=%d, P] = %f\n", s, s_failure, 1.0 - p);
            mdp.setTransitionProbability(s, PLAY, s_failure, 1.0 - p);
        }
    }
    printf ("THIRD CHECK\n");
    mdp.Check();



	
    value_iteration.ComputeStateActionValues(0.01, 10000);
    //value_iteration.ComputeStateValues(0.01, 10000);	

    if (value_iteration.getValue(0, PLAY) > value_iteration.getValue(0, STOP)) {
        printf ("1\n");
    } else {
        printf ("0\n");
    }

    
#if 0
    for (int t=0; t<T; t++) {
        for (int f=0; f<=t; f++) {
            int s = t - f;
            int n = get_bionomial_index(f, s);
            printf ("%+.2f     ", value_iteration.getValue(n));
        }
        printf("\n");
    }
#endif
    for (int t=0; t<T; t++) {
        for (int f=0; f<=t; f++) {
            int s = t - f;
            int n = get_bionomial_index(f, s);
            printf ("%+.4f %+.4f",
                    value_iteration.getValue(n, PLAY),
                    value_iteration.getValue(n, STOP));
        }
        printf("\n");
    }
    return 0;
}

real evaluate_belief_mdp(int T, real gamma, int n_trials, int horizon)
{
    real prior_belief = 0.5;
    real low_value = 0.25;
    real high_value = 0.75;
    real prior = 0.5; // the actual prior


    Vector stats(horizon);
	double init_time = GetCPU();
	double total_time = 0.0;
    for (int trial=0; trial<n_trials; trial++) {
        real actual_value;
        if (urandom()<prior) {
            actual_value = low_value;
        } else {
            actual_value = high_value;
        }
        SimpleBelief initial_belief(low_value, high_value, prior_belief, 1.0 - prior_belief);
        for (int t=0; t<horizon; t++) {
            int a=play_belief_mdp(T, gamma, initial_belief);
            real r = 0.0;
            if (a==0) {
                r = 0.5;
            } else if (a==1) {
                if (urandom() < actual_value) {
                    r = 1.0;
                } else {
                    r = 0.0;
                }
                initial_belief.Observe(r);
            }
            stats[t] += r;
        }
		total_time = GetCPU() - init_time;
		if (trial % 10 == 0) {
			double s_per_trial = total_time / (double) (1+trial);
			
			fprintf(stderr, "T:%.1fs %.1fs/trial ETA:%.1fm\n",
					total_time,
					s_per_trial,
					((double) n_trials - trial)*s_per_trial/60.0);
		}
    }
    real total_reward = 0.0;
    for (int t=0; t<horizon; t++) {
        stats[t] /= (real) n_trials;
        total_reward += stats[t];
        printf("%f\n", stats[t]);
    }
    return total_reward;
}

real test_greedy(int n_trials, int horizon, real alpha)
{
    real low_value = 0.25;
    real high_value = 0.75;
    real prior = 0.5; // the actual prior


    Vector stats(horizon);

    for (int trial=0; trial<n_trials; trial++) {
        real actual_value;
        if (urandom()<prior) {
            actual_value = low_value;
        } else {
            actual_value = high_value;
        }
        real Er[2] = {1.0, 1.0};
        for (int t=0; t<horizon; t++) {
            int a=ArgMax(2, Er);
            real r = 0.0;
            if (a==0) {
                r = 0.5;
            } else if (a==1) {
                if (urandom() < actual_value) {
                    r = 1.0;
                } else {
                    r = 0.0;
                }
            }
            Er[a] += alpha * (r-Er[a]);
            stats[t] += r;
        }
    }
    real total_reward = 0.0;
    for (int t=0; t<horizon; t++) {
        stats[t] /= (real) n_trials;
        total_reward += stats[t];
        printf("%f\n", stats[t]);
    }
    return total_reward;
}


real test_hoeffding_bound(int n_trials, int horizon, real delta)
{
    real low_value = 0.25;
    real high_value = 0.75;
    real prior = 0.5; // the actual prior

    Vector stats(horizon);
    
    int int_max = std::numeric_limits<int>::max();

    real C = 2.0*log(1.0/delta);
    for (int trial=0; trial<n_trials; trial++) {
        real actual_value;
        if (urandom()<prior) {
            actual_value = low_value;
        } else {
            actual_value = high_value;
        }

        int n_actions=2;
        real Er[2] = {0.5, 0.5};
        int N[2] = {int_max, 0};
        real Sr[2] = {0.0, 0.0};

        for (int t=0; t<horizon; t++) {
            for (int i=0; i<n_actions; i++) {
                if (N[i]==int_max) {
                    Sr[i] = 0.0;
                } else if (N[i]==0) {
                    Sr[i] = 1.0;
                } else {
                    Sr[i] = sqrt(C/(real) N[i]);
                }
                Sr[i] += Er[i];
            }
            int a=ArgMax(n_actions, Sr);
            real r = 0.0;
            if (a==0) {
                r = 0.5;
            } else if (a==1) {
                if (urandom() < actual_value) {
                    r = 1.0;
                } else {
                    r = 0.0;
                }
            }
            if (N[a]<int_max) {
                N[a]++;
                Er[a] += (r-Er[a])/((real) N[a]);
            }
            
            stats[t] += r;
        }
    }
    real total_reward = 0.0;
    for (int t=0; t<horizon; t++) {
        stats[t] /= (real) n_trials;
        total_reward += stats[t];
        printf("%f\n", stats[t]);
    }
    return total_reward;
}


int play_belief_mdp(int T, real gamma, SimpleBelief& initial_belief)
{
    int n_states = T*(T+1)/2;
    int n_actions = 2;

    enum Moves {STOP = 0, PLAY};

    real stop_value = 0.5;
    SingularDistribution stop_reward_dist(stop_value);
    SingularDistribution play_one_dist(1.0);
    SingularDistribution play_zero_dist(0.0);


    DiscreteMDP mdp(n_states, n_actions, NULL, NULL);

    ValueIteration value_iteration(&mdp, gamma, 0.5);

    std::vector<SimpleBelief*> belief(n_states);
    std::vector<Pair> pairs(n_states);
    for (int f=0; f<=T; f++) {
        for (int s=0; s<=T; s++) {
            int  t = f+s;
            if (t<T) {
                int n = get_bionomial_index(f, s);
                SMART_ASSERT(n>=0 && n<n_states)(n)(n_states);
                pairs[n].f = f;
                pairs[n].s = s;
                belief[n] = new SimpleBelief(initial_belief);
                for (int i=0; i<f; i++) {
                    belief[n]->Observe(0.0);
                }
                for (int i=0; i<s; i++) {
                    belief[n]->Observe(1.0);
                }
            }
        }
    }

    // First zero all transitions
    for (int s=0; s<n_states; s++) {
        for (int s2=0; s2<n_states; s2++) {
            mdp.setTransitionProbability(s, PLAY, s2, 0.0);
            mdp.setTransitionProbability(s, STOP, s2, 0.0);
        }
    }
    //printf ("FIRST CHECK\n");
    //mdp.Check();
	
    // When stopping, you don't get to observe nothin'
    for (int s=0; s<n_states; s++) {
        mdp.setRewardDistribution(s, STOP, &stop_reward_dist);
        mdp.setTransitionProbability(s, STOP, s, 1.0);
    }

    //printf ("SECOND CHECK\n");
    //mdp.Check();
    // When playing, beliefs change, until the end of time, where they are
    // stationary.
    std::vector<Distribution*> tmp_dist;
    for (int s=0; s<n_states; s++) {
        float p = belief[s]->getMean();
        Pair b=pairs[s];
        int s_success = get_bionomial_index(b.f, b.s + 1);
        int s_failure = get_bionomial_index(b.f + 1, b.s);
        BernoulliDistribution* dist = new BernoulliDistribution(p);
        tmp_dist.push_back(dist);
        mdp.setRewardDistribution(s, PLAY, dist);
        if (s_success >= n_states || s_failure >= n_states) {
            mdp.setTransitionProbability(s, PLAY, s, 1.0);
        } else {
            mdp.setTransitionProbability(s, PLAY, s_success, p);
            mdp.setTransitionProbability(s, PLAY, s_failure, 1.0 - p);
        }
    }
    //printf ("THIRD CHECK\n");
    //mdp.Check();

	
	
    value_iteration.ComputeStateActionValues(0.1, 100);
    //value_iteration.ComputeStateValues(0.01, 10000);	

    for (int n=0; n<n_states; n++) {
        delete belief[n];
    }
    for (unsigned int n=0; n<tmp_dist.size(); n++) {
        delete tmp_dist[n];
    }
    if (value_iteration.getValue(0, PLAY) > value_iteration.getValue(0, STOP)) {
        return 1;
    } else {
        return 0;
    }
}

int play_belief_mdp_quick(int T, real gamma, SimpleBelief& initial_belief)
{
    int n_states = T*(T+1)/2 + 1;
    int n_actions = 2;
    int terminal_state = n_states - 1;

    enum Moves {STOP = 0, PLAY};

    real stop_value = 0.5;

    // Maybe we should take T into account here?
    SingularDistribution stop_reward_dist(stop_value / (1.0 - gamma));
    SingularDistribution play_one_dist(1.0);
    SingularDistribution play_zero_dist(0.0);
    

    DiscreteMDP mdp(n_states, n_actions, NULL, NULL);

    ValueIteration value_iteration(&mdp, gamma, 0.5);

    std::vector<SimpleBelief*> belief(n_states);
    std::vector<Pair> pairs(n_states);
    for (int f=0; f<=T; f++) {
        for (int s=0; s<=T; s++) {
            int  t = f+s;
            if (t<T) {
                int n = get_bionomial_index(f, s);
                SMART_ASSERT(n>=0 && n<n_states)(n)(n_states);
                pairs[n].f = f;
                pairs[n].s = s;
                belief[n] = new SimpleBelief(initial_belief);
                for (int i=0; i<f; i++) {
                    belief[n]->Observe(0.0);
                }
                for (int i=0; i<s; i++) {
                    belief[n]->Observe(1.0);
                }
            }
        }
    }

    // First zero all transitions
    for (int s=0; s<n_states; s++) {
        for (int s2=0; s2<n_states; s2++) {
            mdp.setTransitionProbability(s, PLAY, s2, 0.0);
            mdp.setTransitionProbability(s, STOP, s2, 0.0);
        }
    }
    //printf ("FIRST CHECK\n");
    //mdp.Check();
	
    // When stopping, you don't get to observe nothin'
    for (int s=0; s<n_states; s++) {
        mdp.setRewardDistribution(s, STOP, &stop_reward_dist);
        mdp.setTransitionProbability(s, STOP, terminal_state, 1.0);
    }

    mdp.setTransitionProbability(terminal_state, PLAY, terminal_state, 1.0);
    mdp.setRewardDistribution(terminal_state, PLAY, &stop_reward_dist);
    mdp.setRewardDistribution(terminal_state, STOP, &stop_reward_dist);

    //printf ("SECOND CHECK\n");
    //mdp.Check();
    // When playing, beliefs change, until the end of time, where they are
    // stationary.
    // Set distribution for non-terminal states
    std::vector<Distribution*> tmp_dist;
    for (int s=0; s<terminal_state; s++) {
        float p = belief[s]->getMean();
        Pair b=pairs[s];
        int s_success = get_bionomial_index(b.f, b.s + 1);
        int s_failure = get_bionomial_index(b.f + 1, b.s);
        BernoulliDistribution* dist = new BernoulliDistribution(p);
        tmp_dist.push_back(dist);
        mdp.setRewardDistribution(s, PLAY, dist);
        if (s_success >= n_states || s_failure >= n_states) {
            mdp.setTransitionProbability(s, PLAY, s, 1.0);
        } else {
            mdp.setTransitionProbability(s, PLAY, s_success, p);
            mdp.setTransitionProbability(s, PLAY, s_failure, 1.0 - p);
        }
    }
    //printf ("THIRD CHECK\n");
    //mdp.Check();

	
	
    value_iteration.ComputeStateActionValues(0.01, 10000);
    //value_iteration.ComputeStateValues(0.01, 10000);	

    for (int n=0; n<n_states; n++) {
        delete belief[n];
    }
    for (unsigned int n=0; n<tmp_dist.size(); n++) {
        delete tmp_dist[n];
    }
    if (value_iteration.getValue(0, PLAY) > value_iteration.getValue(0, STOP)) {
        return 1;
    } else {
        return 0;
    }
}



#endif

