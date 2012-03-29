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


class ObjectiveBanditPolicy
{
public:
	const int n_actions;
	const int n_outcomes;
	ObjectiveBanditPolicy(int n_actions_, int n_outcomes_)
		: n_actions(n_actions_), n_outcomes(n_outcomes_)
	{
	}
	virtual ~ObjectiveBanditPolicy()
	{}
	virtual int getAction(const Vector& payoff) const = 0;
	virtual int Act(const Vector& payoff, int outcome) = 0;
};



class GreedyObjectiveBandit : public ObjectiveBanditPolicy
{
protected:
	int action;
public:
	Matrix N; ///< probability matrix
	Matrix P; ///< probability matrix
	GreedyObjectiveBandit(int n_actions_, int n_outcomes_)
		: ObjectiveBanditPolicy(n_actions_, n_outcomes_),
		  action(-1),
		  N(n_actions, n_outcomes),
		  P(n_actions, n_outcomes)
	{
		real p = 1.0 / (real) n_outcomes;
		for (int i=0; i<n_actions; ++i) {
			for (int j=0; j<n_outcomes; ++j) {
				N(i,j) = 0.5;
				P(i,j) = p;
			}
		}
	}

	virtual ~GreedyObjectiveBandit()
	{
	}

	virtual int getAction(const Vector& payoff) const
	{
		const Matrix& rP = P;
		return ArgMax(rP * payoff);
	}
	virtual int Act(const Vector& payoff, int outcome)
	{
		if (action >= 0 && outcome >= 0) {
			N(action, outcome)++;
			Vector p = N.getRow(action);
			P.setRow(action, p / p.Sum());
		}
		action = getAction(payoff);

		return action;
	}
};

struct RegretPair
{
	real	average;
	real worst;
};

RegretPair getRegret(const Matrix& P, const GreedyObjectiveBandit& policy)
{
	RegretPair regret;
	regret.average = 0.0;
	regret.worst = 0.0;
	int n_actions = policy.n_actions;
	int n_outcomes = policy.n_outcomes;

	for (int i=0; i<n_outcomes; ++i) {
		Vector payoff(n_outcomes);
		payoff(i) = 1.0;
		int action = policy.getAction(payoff);
		const Vector& r_payoff = payoff;
		Vector U = P * r_payoff;
		real delta = Max(U) - U(action);
		regret.average += delta;
		regret.worst = std::max(regret.worst, delta);
	}
	regret.average /= (real) n_outcomes;
	return regret;
}




int main (int argc, char** argv)
{

    srand48(1228517343);
    //RandomNumberFile rng("./dat/r1e7.bin");
	int n_actions = 4;
	int n_outcomes = 4;
	int horizon = 1000;
	int n_experiments = 10;

	if (argc <= 4) {
		fprintf(stderr, "arguments: n_actions n_outcomes horizon n_runs\n");
		exit(-1);
	}
	n_actions = atoi(argv[1]);
	n_outcomes = atoi(argv[2]);
	horizon = atoi(argv[3]);
	n_experiments = atoi(argv[4]);

	Vector worst_case(horizon);
	Vector average_case(horizon);
    for (int experiment=0; experiment<n_experiments; experiment++) {

		Matrix P_sa(n_actions, n_outcomes);
		DirichletDistribution prior(n_outcomes);
		for (int i=0; i<n_actions; ++i) {
			P_sa.setRow(i, prior.generate());
		}
		GreedyObjectiveBandit policy(n_actions, n_outcomes);
		real reward = 0.0;
		int outcome = -1;
        for (int t=0; t<horizon; t++) {
			Vector payoff = prior.generate();
            int action = policy.Act(payoff, outcome);
			outcome = DiscreteDistribution::generate(P_sa.getRow(action));
			reward = payoff(outcome);
			RegretPair regret =  getRegret(P_sa, policy);
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
