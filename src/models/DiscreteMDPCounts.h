
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
#include "Dirichlet.h"
#include "DirichletFiniteOutcomes.h"
#include "MeanEstimator.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"
#include "real.h"
#include <vector>

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
    std::vector<DirichletDistribution> P; ///< Vector of Dirichlet distributions on P
    //std::vector<DirichletFiniteOutcomes> P; ///< Vector of Dirichlet distributions on P
    std::vector<ConjugatePrior*> ER; ///< Vector of estimators on ER.
    //std::vector<BetaDistribution> ER; ///< Vector of estimators on ER.
    //std::vector<NormalUnknownMeanPrecision> ER; ///< Vector of estimators on ER.
    DiscreteMDP mean_mdp; ///< a model of the mean MDP
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
        int ID = getID(s, a);
        return (int) ceil(P[ID].getMass());
    }
    //void SetNextReward(int s, int a, real r);
};


#endif

