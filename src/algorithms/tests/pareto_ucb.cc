// -*- Mode: C++; -*-
// copyright (c) 2012 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "Matrix.h"
#include "Dirichlet.h"
#include "Sarsa.h"
#include "MersenneTwister.h"
#include "Bounds.h"
#include "ObjectiveGeneralBandit.h"
#include "ProductDistribution.h"
#include "GammaDistribution.h"
#include "BetaDistribution.h"
#include <algorithm>

struct RegretPair
{
	real	average;
	real worst;
};

RegretPair getRegret(const Matrix& P, const ObjectiveGeneralBandit* policy)
{
	RegretPair regret;
	regret.average = 0.0;
	regret.worst = 0.0;
	//int n_actions = policy.n_actions;
	int n_outcomes = policy->n_outcomes;

	for (int i=0; i<n_outcomes; ++i) {
		Vector payoff(n_outcomes);
		payoff(i) = 1.0;
		int action = policy->getGreedyAction(payoff);
		const Vector& r_payoff = payoff;
		Vector U = P * r_payoff;
		real delta = Max(U) - U(action);
		regret.average += delta;
		regret.worst = std::max(regret.worst, delta);
	}
	regret.average /= (real) n_outcomes;
	return regret;
}


enum Method {
	UNDEFINED = 0x0,
	EPSILON_GREEDY,
    HOEFFDING_UCB
};

int main (int argc, char** argv)
{

    srand48(1228517343);
	MersenneTwisterRNG rng;
	rng.manualSeed(1228517343);
	setRandomSeed(1228517343);
    //RandomNumberFile rng("./dat/r1e7.bin");
	int n_actions = 4;
	int n_outcomes = 4;
	int horizon = 1000;
	int n_experiments = 10;
	real epsilon = 0.0;
	Method method = UNDEFINED;

	if (argc <= 6) {
		fprintf(stderr, "arguments: n_actions n_outcomes horizon n_runs method epsilon \n method 1: epsilon-greedy \n method 2: hoeffding_ucb\n");
		exit(-1);
	}
	n_actions = atoi(argv[1]);
	n_outcomes = atoi(argv[2]);
	horizon = atoi(argv[3]);
	n_experiments = atoi(argv[4]);
	method = (Method) atoi(argv[5]);
	epsilon = atof(argv[6]);
	
	Vector worst_case(horizon);
	Vector average_case(horizon);



    BetaDistribution outcome_parameter_prior(2, 2);	

    for (int experiment=0; experiment<n_experiments; experiment++) {
        
		Matrix P_x_a(n_actions, n_outcomes);
		DirichletDistribution payoff_prior(n_outcomes);

		for (int i=0; i<n_actions; ++i) {
            for (int j=0; j<n_outcomes; ++j) {
                P_x_a(i, j) = outcome_parameter_prior.generate();
            }
		}

		ObjectiveGeneralBandit* policy = NULL;
		switch (method) {
		case EPSILON_GREEDY:
			policy = new EpsilonGreedyObjectiveGeneralBandit(n_actions, n_outcomes, rng, epsilon);
			break;
		case HOEFFDING_UCB:
			policy = new HoeffdingObjectiveGeneralBandit(n_actions, n_outcomes, rng);
			break;
		default:
			fprintf(stderr, "Unknown method %d\n", method);
			exit(-1);
		}	
		real reward = 0.0;
        Vector outcome(n_outcomes);
        for (int t=0; t<horizon; t++) {
			Vector payoff = payoff_prior.generate();
            int action = policy->Act(payoff, outcome);
            for (int j=0; j<n_outcomes; ++j) {
                if (urandom() < P_x_a(action, j)) {
                    outcome(j) = 1.0;
                } else {
                    outcome(j) = 0.0;
                }
            }
			reward = Product(payoff, outcome);
			RegretPair regret =  getRegret(P_x_a, policy);
			worst_case(t) += regret.worst;
			average_case(t) += regret.average;
		}
    }

	worst_case /= (real) n_experiments;
	average_case /= (real) n_experiments;
	for (int t=0; t<horizon; ++t) {
		printf ("%f %f\n", worst_case(t), average_case(t));
	}

    return 0;
}



#endif
