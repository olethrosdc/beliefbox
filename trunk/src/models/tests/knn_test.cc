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

bool knn_test(int n_neighbours, real rbf_beta)
{
    int n_actions, n_dim;
    real T = 100.0;
    real dt = 0.01;
    n_actions = 1;
    n_dim = 2;

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
        TrajectorySample sample(s, 0, 0.0, s2);
        model.AddSample(sample);
        model.GetExpectedTransition(s, 0, predicted_reward, predicted_state, n_neighbours, rbf_beta);
        real err = EuclideanNorm(&s2, &predicted_state);
        printf ("%f %f %f\n", s[0], s[1], err);
        s = s2;
        model.ValueIteration(n_neighbours, rbf_beta);
    }
    return true;
}

bool knn_environment_test(Environment<Vector, int>& environment, int n_neighbours, real rbf_beta, int T)
{
    int n_actions, n_dim;
    n_actions = 3;
    n_dim = 2;


    KNNModel model(n_actions, n_dim, 0.99);
    
    Vector state(n_dim);
    Vector next_state(n_dim);
    real reward = 0.0;
    Vector predicted_state(n_dim);
    //real predicted_reward = 0.0;

    environment.Reset();

    int action = rand()%n_actions;
    for (int t=0; t<T; ++t) {
        state = environment.getState();
        Vector state_copy = state;
        action = model.GetBestAction(state_copy, n_neighbours, rbf_beta);
        bool action_ok = environment.Act(action);


        next_state = environment.getState();
        reward = environment.getReward();
        TrajectorySample sample(state, action, reward, next_state);
        model.AddSample(sample);
        state_copy = state;
        //        model.GetExpectedTransition(state_copy, action, predicted_reward, predicted_state, n_neighbours, rbf_beta);        
        //real err = EuclideanNorm(&next_state, &predicted_state);
        //printf ("%f %f %f %f %d\n", state[0], state[1], err, reward, action);
        real V = model.GetExpectedValue(state_copy, n_neighbours, rbf_beta);
        printf ("%f %f %f\n", state[0], state[1], V);
        if (action_ok) {
            model.ValueIteration(n_neighbours, rbf_beta);
        } else {
            printf ("# BEEP\n");
            environment.Reset();
        }
    }
    
    for (real x = -1; x<1; x += 0.01) {
        for (real y = -1; y<1; y += 0.01) {
            state[0] = x;
            state[1] = y;
            real V = model.GetExpectedValue(state, n_neighbours, rbf_beta);
            printf ("%f ",  V);
        }
        printf ("# V\n");
    }
    return true;
}

int main (int argc, char** argv)
{
    printf ("# knn_test.cc: testing model adaptation\n");
    if (argc != 3) {
        fprintf (stderr, "usage: knn_test n_neighbours rbf_beta\n");
        exit(-1);
    }
    bool success = true;//knn_test(5, 1.0);
    int n_neighbours = atoi(argv[1]);
    real rbf_beta = atof(argv[2]);
    printf ("# neighbours: %d - RBF beta: %f\n", n_neighbours, rbf_beta);

    MountainCar mountain_car;
    real start_time = GetCPU();
    success = knn_environment_test(mountain_car, n_neighbours, rbf_beta, 100);
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
