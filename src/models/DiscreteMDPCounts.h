
// -*- Mode: c++ -*-
// copyright (c) 2005-2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DISCRETE_MDP_COUNTS
#define DISCRETE_MDP_COUNTS

#include "MDPModel.h"
#include "DirichletTransitions.h"
#include "MeanEstimator.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"
#include "real.h"
#include <vector>
#include <unordered_map>

/** This implementation of an MDP model is based on transition counts.
 */
class DiscreteMDPCounts : public MDPModel
{
public:
    enum RewardFamily {
        UNDEFINED=0x0,
        BETA,
        NORMAL,
        FIXED
    };
protected:
	bool use_sampling = false;
	/// Dirichlet distribution for transitions
    DirichletTransitions transitions; 
	/// Vector of estimators on ER.
    std::vector<ConjugatePrior*> ER; 
    DiscreteMDP mean_mdp; ///< a model of the mean MDP
	DiscreteMDP* sampled_mdp = NULL; ///< a model of the mean MDP 
    RewardFamily reward_family; ///< reward family to be used
    int N;
    int getID (int s, int a) const
    {
        assert(s>=0 && s<n_states);
        assert(a>=0 && a<n_actions);
        return s*n_actions + a;
    }
    Vector getDirichletParameters (int s, int a) const;
public:
    DiscreteMDPCounts (int n_states, int n_actions, real init_transition_count= 0.5, RewardFamily reward_family=NORMAL);
	void useSampling(bool sampling) {
		use_sampling = sampling;
		if (!sampled_mdp) {
			sampled_mdp = generate();
		}
	}
	DiscreteMDPCounts(const DiscreteMDPCounts& model) :
		use_sampling(model.use_sampling),
		transitions(model.transitions),
		mean_mdp(model.mean_mdp),
		reward_family(model.reward_family)
		N(model.N)
	{
		//printf("Creating DiscreteMDPCounts with %d states and %d actions\n",  n_states, n_actions);
		ER.resize(N);
		for (int i=0; i<N; ++i) {
			switch(reward_family) {
			case BETA:
				ER[i] = new BetaDistribution();
				break;
			case NORMAL:
				ER[i] = new NormalUnknownMeanPrecision();
				break;
			case FIXED:
				ER[i] = new UnknownSingularDistribution();
				break;
			default:
				Serror("Unknown distribution family %d\n", reward_family);
			}
		}
		for (int s=0; s<n_states; s++) {
			for (int a=0; a<n_actions; a++) {
				for (int s_next=0; s_next<n_states; s_next++) {
					real p = transitions.marginal_pdf(s, a, s_next);
					mean_mdp.setTransitionProbability(s, a, s_next, p);
					real expected_reward = getExpectedReward(s,a);
					mean_mdp.reward_distribution.setFixedReward(s, a, expected_reward);
				}
			}
		}
	}
	virtual DiscreteMDPCounts* Clone();
    virtual ~DiscreteMDPCounts();
    virtual void AddTransition(int s, int a, real r, int s2);
    virtual void setFixedRewards(const Matrix& rewards);
    virtual real GenerateReward (int s, int a) const;
    virtual int GenerateTransition (int s, int a) const;
    virtual real getTransitionProbability (int s, int a, int s2) const;
    virtual Vector getTransitionProbabilities (int s, int a) const;
    virtual real getExpectedReward (int s, int a) const;

    virtual void Reset();
    virtual void ShowModel() const;

    virtual DiscreteMDP* generate() const;
    virtual const DiscreteMDP* const getMeanMDP() const;
    //virtual DiscreteMDP* CreateMDP() const;
    virtual void CopyMeanMDP(DiscreteMDP* mdp) const;
    int getNVisits(int s, int a) const
    {
		return transitions.getCounts(s, a);
    }
    //void SetNextReward(int s, int a, real r);
};


#endif

