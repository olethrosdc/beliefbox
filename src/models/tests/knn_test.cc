/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "KNNModel.h"
#include "Vector.h"
#include "NormalDistribution.h"
#include "Random.h"
#include "Vector.h"
#include "EasyClock.h"
#include "MountainCar.h"
#include "Pendulum.h"
#include "ContinuousChain.h"

bool knn_test(int n_neighbours, real rbf_beta)
{
    int n_actions, n_dim;
    real T = 100.0;
    real dt = 0.01;
    n_actions = 1;
    n_dim = 2;

    real alpha = 0.5;

    KNNModel model(n_actions, n_dim);
    
    Vector s(n_dim);
    Vector s2(n_dim);
    Vector predicted_state(n_dim);
    real predicted_reward;
    real u = 1;
    real x =  0.0;
    real K = 10.0;
    real b = 0.0001;
    real nu = 0.01;
    s[0] = x;
    s[1] = u;
    for (real t=0; t<=T; t += dt) {
        real a = - b * u - K * x + urandom(0, nu);
        u += a * dt;
        x += u * dt;
        s2[0] = x;
        s2[1] = u;
        TrajectorySample sample(s, 0, 0.0, s2, false);
        model.AddSample(sample);
        model.GetExpectedTransition(alpha, s, 0, predicted_reward, predicted_state, n_neighbours, rbf_beta);
        real err = EuclideanNorm(&s2, &predicted_state);
        printf ("%f %f %f\n", s[0], s[1], err);
        s = s2;
        model.ValueIteration(alpha, n_neighbours, rbf_beta);
    }
    return true;
}

bool knn_environment_test(Environment<Vector, int>& environment, real alpha, int n_neighbours, real rbf_beta, real gamma, int T, int max_samples)
{
    int n_actions = environment.getNActions();
    int n_dim = environment.getNStates();
    printf ("# Starting test with %d state dimensions and %d actions.\n", n_dim, n_actions);

    KNNModel model(n_actions, n_dim, gamma, false, 0.0001, 1.0);
    model.SetMaxSamples(max_samples);

    Vector state(n_dim);
    Vector next_state(n_dim);
    real reward = 0.0;
    Vector predicted_state(n_dim);
    real predicted_reward = 0.0;

    environment.Reset();

    int action = rand()%n_actions;
    int n_steps = 0;
    for (int t=0; t<T; ++t) {
        n_steps++;
        // select an action
        real epsilon = 1000.0 / (1000.0 + (real) t);

        state = environment.getState();
        if (t < n_actions) {
            action = t;
        } else {
            if (urandom() < epsilon) {
                action = rand()%n_actions;
            } else {
                action = model.GetBestAction(state, n_neighbours, rbf_beta);
            }
        }


        // act
        bool action_ok = environment.Act(action);
        next_state = environment.getState();
        reward = environment.getReward();

        // update model
        TrajectorySample sample(state, action, reward, next_state, !action_ok);
        model.AddSample(sample);
        model.UpdateValue(sample, alpha, n_neighbours, rbf_beta);
        // predict next state
        model.GetExpectedTransition(alpha, state, action, predicted_reward, predicted_state, n_neighbours, rbf_beta);        
        real err = EuclideanNorm(&next_state, &predicted_state);
        for (int i=0; i<=10; ++i) {
            model.ValueIteration(alpha, n_neighbours, rbf_beta);
        }
        // update value function
        real V = model.GetExpectedValue(state, n_neighbours, rbf_beta);
        if (n_dim == 2) {
            printf ("%f %f %d %f %f %f %f\n", state[0], state[1], action, epsilon, V, err, reward);
        } else {
            printf ("%f %d %f %f %f %f\n", state[0], action, epsilon, V, err, reward);
        }

        // see if we are successful
        if (!action_ok) {
            for (int a=0; a<=n_actions; ++a) {
                TrajectorySample sample(environment.getState(), action, reward, environment.getState(), true);
            }
            model.AddSample(sample);

            //printf ("%d # steps: %f %f\n", n_steps, state[0], state[1]);

            n_steps = 0;
            
            environment.Reset();
        }
    }
    
    printf ("## final round of value iteration\n");
    for (int i=0; i<10000; i++) {
        model.ValueIteration(alpha, n_neighbours, rbf_beta);
    }    

    if (n_dim == 1) {
        real x_min = -1.5;
        real x_max = 1.5;
        real dx = (x_max - x_min) / 100.0;
        for (real x = x_min; x<x_max; x += dx) {
            state[0] = x;
            real V = model.GetExpectedValue(state, n_neighbours, rbf_beta);
            printf ("%f %f # V\n",  x, V);
        }
        printf ("# V\n");

        for (real x = x_min; x<x_max; x += dx) {
            state[0] = x;
            printf ("%f ", x);
            for (int a=0; a<n_actions; ++a) {
                real V = model.GetExpectedActionValue(state, a, n_neighbours, rbf_beta);
                printf ("%f ",  V);
            }
            printf ("# Q\n");
        }

        for (real x = x_min; x<x_max; x += dx) {
            Vector next_state(n_dim);
            real next_reward;
            state[0] = x;
            printf ("%f ", x);
            for (int a=0; a<n_actions; ++a) {
                //model.GetExpectedTransition(alpha, state, action, predicted_reward, predicted_state, n_neighbours, rbf_beta);        
                model.GetExpectedTransition(alpha, state, a, next_reward, next_state,n_neighbours, rbf_beta);
                printf ("%f ",  next_state[0]);
            }
            printf ("# STATE\n");
        }
    } else if (n_dim == 2) {
    // print out value function
        real x_min = -1.5;
        real x_max = 1.5;
        real y_min = -3.0;
        real y_max = 3.0;
        real dx = (x_max - x_min) / 100.0;
        real dy = (y_max - y_min) / 100.0;
        
        for (real x = x_min; x<x_max; x += dx) {
            for (real y = y_min; y<y_max; y += dy) {
                state[0] = x;
                state[1] = y;
                real V = model.GetExpectedValue(state, n_neighbours, rbf_beta);
                printf ("%f ",  V);
            }
            printf ("# V\n");
        }
        
        for (real x = x_min; x<x_max; x += dx) {
            for (real y = y_min; y<y_max; y += dy) {
                state[0] = x;
                state[1] = y;
                for (int a=0; a<n_actions; ++a) {
                    real V = model.GetExpectedActionValue(state, a, n_neighbours, rbf_beta);
                    printf ("%f ",  V);
                }
                printf ("# Q\n");
            }
        }
    }

    

    model.Show();

    return true;
}

int main (int argc, char** argv)
{
    printf ("# knn_test.cc: testing model adaptation\n");
    if (argc != 8) {
        fprintf (stderr, "usage: knn_test task alpha n_neighbours rbf_beta gamma horizon n_samples\n");
        fprintf (stderr, "task in {chain, pendulum, car}\n alpha in [0,1]\n n_neighbours in {1, 2, ...}\n rbf_beta > 0\n");
        exit(-1);
    }
    bool success = true;//knn_test(5, 1.0);
    char* task = argv[1];
    real alpha = atof(argv[2]);
    int n_neighbours = atoi(argv[3]);
    real rbf_beta = atof(argv[4]);
    real gamma = atof(argv[5]);
    int horizon = atoi(argv[6]);
    int n_samples = atoi(argv[7]);
    printf ("# neighbours: %d - RBF beta: %f\n", n_neighbours, rbf_beta);

    MountainCar mountain_car;
    Pendulum pendulum;
    ContinuousChain continuous_chain;
    real start_time = GetCPU();
    if (!strcmp(task, "car")) {
        success = knn_environment_test(mountain_car, alpha, n_neighbours, rbf_beta, gamma, horizon, n_samples);
    } else  if (!strcmp(task, "pendulum")) {
        success = knn_environment_test(pendulum, alpha, n_neighbours, rbf_beta, gamma, horizon, n_samples);
    } else if (!strcmp(task, "chain")) {
        success = knn_environment_test(continuous_chain, alpha, n_neighbours, rbf_beta, gamma, horizon, n_samples);
    }
    real end_time = GetCPU();
    printf ("# Time: %f\n", end_time - start_time);

    if (success)  {
        printf ("# Test succeeded\n");
        return 0;
    }
    
    fprintf (stderr, "# Test failed\n");
    return -1;
}

#endif
