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

#include "BasisSet.h"

RBFBasisSet::RBFBasisSet(const EvenGrid& grid, real scale)
{
    n_bases = 0;
    for (int i=0; i<grid.getNIntervals(); ++i) {
        AddCenter(grid.getCenter(i), grid.delta * scale);
    }
    //logmsg("Added %d RBFs\n", centers.size());
}

RBFBasisSet::~RBFBasisSet()
{
    for (uint i=0;  i<centers.size(); ++i) {
        delete centers[i];
    }
}

void RBFBasisSet::AddCenter(const Vector& v, const Vector& b)
{
    RBF* rbf = new RBF(v, b);
    centers.push_back(rbf);
    n_bases++;
    features.Resize(centers.size());
    features[n_bases-1] = 0.0;
    log_features.Resize(centers.size());
    log_features[n_bases-1] = 0.0;
    valid_features = false;
    valid_log_features = false;
}

void RBFBasisSet::AddCenter(const Vector& v, real b)
{
    RBF* rbf = new RBF(v, b);
    centers.push_back(rbf);
    n_bases++;
    features.Resize(centers.size());
    features[n_bases-1] = 0.0;
    log_features.Resize(centers.size());
    log_features[n_bases-1] = 0.0;
    valid_features = false;
    valid_log_features = false;
}

void RBFBasisSet::logEvaluate(const Vector& x) const
{
    real log_sum = LOG_ZERO;
    for (int i=0; i<n_bases; ++i) {
        log_features[i] = centers[i]->logEvaluate(x);
        log_sum = logAdd(log_features[i], log_sum);
    }
    for (int i=0; i<n_bases; ++i) {
        log_features[i] -= log_sum;
    }
    valid_log_features = true;
    valid_features = false;
}


void RBFBasisSet::Evaluate(const Vector& x) const
{
    for(int i = 0; i<n_bases; ++i){
        features[i] = centers[i]->Evaluate(x);
    }
    valid_log_features = false;
    valid_features = true;
}

void RBFBasisSet::Observe(const int& action, real reward, const Vector& next_state)
{
    logEvaluate(next_state);
    for (int i=0; i<n_bases; ++i) {
        features[i] = exp(log_features[i]);
    }
	valid_log_features = true;
    valid_features = true;
}

void RBFBasisSet::Observe(const Vector& state)
{
    logEvaluate(state);
    for (int i=0; i<n_bases; ++i) {
        features[i] = exp(log_features[i]);
    }
	valid_log_features = true;
    valid_features = true;
}

void CountsBasis::Reset()
{
	state = -1;
}

CountsBasis::~CountsBasis()
{
	// nothing to do
}


void CountsBasis::Observe(const int& action, real reward, const int& next_state)
{
	if (state >= 0) {
		int N = model.getCounts(state, action);
		real alpha = 1.0f / (real) (N + 1);
		average_reward(state, action) += alpha * (reward - 			average_reward(state, action));
	}
			
	model.Observe(state, action, next_state);
	
	// split observation in three parts.
	// 1: current state |S|
	// 2: current reward 1
	// 3: average transitions |S^2A|
	// 4: transition counts |S^2A|
	// 5: empirical average rewards |SA|
	features.Clear();
	features(next_state) = 1;
	features(n_states) = reward;
	int i = n_states + 1;


	
	for (int s = 0; s < n_states; ++s) {
		for (int a = 0; a < n_actions; ++a) {
			Vector marginal = model.getMarginal(s, a);
			for (int s2 = 0; s2 < n_states; ++s2, ++i) {
				features(i) = marginal(s2);
			}
		}
	}

	for (int s = 0; s < n_states; ++s) {
		for (int a = 0; a < n_actions; ++a) {
			Vector counts = model.getParameters(s, a);
			for (int s2 = 0; s2 < n_states; ++s2, ++i) {
				features(i) = log(1 + counts(s2));
			}
		}
	}


	for (int s = 0; s < n_states; ++s) {
		for (int a = 0; a < n_actions; ++a, ++i) {
			features(i) = average_reward(s, a);
		}
	}
	state = next_state;
	valid_features = true;
	valid_log_features = false;
}
