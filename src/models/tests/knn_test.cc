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
    }
    return true;
}

bool knn_mountain_car_test(int n_neighbours, real rbf_beta)
{
    int n_actions, n_dim;
    real T = 100.0;
    real dt = 0.01;
    n_actions = 2;
    n_dim = 2;


    KNNModel model(n_actions, n_dim);
    
    Vector s(n_dim);
    Vector s2(n_dim);
    real r = 0.0;
    Vector predicted_state(n_dim);
    real predicted_reward;

    MountainCar environment;
    environment.Reset();
    s = environment.getState();
    for (real t=0; t<=T; t += dt) {
        int action = rand()%n_actions;
        bool action_ok = environment.Act(action);
        s2 = environment.getState();
        r = environment.getReward();
        TrajectorySample sample(s, action, r, s2);
        model.AddSample(sample);
        model.GetExpectedTransition(s, 0, predicted_reward, predicted_state, n_neighbours, rbf_beta);
        real err = EuclideanNorm(&s2, &predicted_state);
        printf ("%f %f %f\n", s[0], s[1], err);
        s = s2;
    }
    return true;
}

int main (int argc, char** argv)
{
    printf ("# knn_test.cc: testing model adaptation\n");
    bool success = knn_test(5, 1.0);
    success = (success && knn_mountain_car_test(5, 1.0));
    if (success)  {
        printf ("# Test succeeded\n");
        return 0;
    }
    
    fprintf (stderr, "# Test failed\n");
    return -1;
}

#endif
