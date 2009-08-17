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

bool knn_environment_test(Environment<Vector, int>& environment, real alpha, int n_neighbours, real rbf_beta, int T)
{
    int n_actions, n_dim;
    n_actions = 3;
    n_dim = 2;

    real gamma = 1.0;
    KNNModel model(n_actions, n_dim, gamma, true, 0.001, 0.0);
    
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
        real epsilon = 100.0 / (100.0 + (real) t);

        state = environment.getState();
        if (t < n_actions) {
            action = t;
        } else {
            action = model.GetBestAction(state, n_neighbours, rbf_beta, epsilon);
        }
        if (t < 10) {
            action = 1;
            if (state[0] + 0.1* state[1] > 0.125) {
                action = 2;
            } else if (state[0] + 0.1*state[1] < -0.125) {
                action = 0;
            }
        }



        // act
        bool action_ok = environment.Act(action);
        next_state = environment.getState();
        reward = environment.getReward();

        // update model
        TrajectorySample sample(state, action, reward, next_state, !action_ok);
        if (t < 1000) {
            model.AddSample(sample);
        }

        // predict next state
        model.GetExpectedTransition(alpha, state, action, predicted_reward, predicted_state, n_neighbours, rbf_beta);        
        real err = EuclideanNorm(&next_state, &predicted_state);

        // update value function
        real V = model.GetExpectedValue(state, n_neighbours, rbf_beta);
        printf ("%f %f %d %f %f %f %f\n", state[0], state[1], action, epsilon, V, err, reward);

        // see if we are successful
        if (action_ok) {
            for (int iter=0; iter<1; ++iter) {
                model.ValueIteration(alpha, n_neighbours, rbf_beta);
            }
        } else {
            printf ("%d # steps: %f %f\n", n_steps, state[0], state[1]);
            n_steps = 0;
            environment.Reset();
        }
    }
    
    
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

    for (int a=0; a<n_actions; ++a) {
        for (real x = x_min; x<x_max; x += dx) {
            for (real y = y_min; y<y_max; y += dy) {
                state[0] = x;
                state[1] = y;
                real V = model.GetExpectedActionValue(state, a, n_neighbours, rbf_beta);
                printf ("%f ",  V);
            }
            printf ("# Q%d\n", a);
        }
    }
    return true;
}

int main (int argc, char** argv)
{
    printf ("# knn_test.cc: testing model adaptation\n");
    if (argc != 4) {
        fprintf (stderr, "usage: knn_test alpha n_neighbours rbf_beta\n");
        exit(-1);
    }
    bool success = true;//knn_test(5, 1.0);
    real alpha = atof(argv[1]);
    int n_neighbours = atoi(argv[2]);
    real rbf_beta = atof(argv[3]);
    printf ("# neighbours: %d - RBF beta: %f\n", n_neighbours, rbf_beta);

    MountainCar mountain_car;
    Pendulum pendulum;
    real start_time = GetCPU();
    success = knn_environment_test(mountain_car, alpha, n_neighbours, rbf_beta, 10000);
    //success = knn_environment_test(pendulum, alpha, n_neighbours, rbf_beta, 20000);
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
