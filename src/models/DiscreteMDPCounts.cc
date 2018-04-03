// -*- Mode: c++ -*-
// copyright (c) 2005-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "DiscreteMDPCounts.h"
#include "Random.h"
#include "SingularDistribution.h"

#include <stdexcept>

/** Create a counting model of an MDP.
	
	\arg n_states the number of MDP states
	\arg n_actions the number of MDP actions
	\arg init_transition_count the prior for the Dirichlet. The higher this is, the more the model will expect to see unseen transitions.
	\arg init_reward_count the prior number of counts for the reward.  This should be used in conjuction with the prior reward average to bias the rewards.
	\arg The prior reward average. This can be used to bias the average to some particular value.
*/
DiscreteMDPCounts::DiscreteMDPCounts (int n_states, int n_actions, real init_transition_count, RewardFamily reward_family_)
    : 
    MDPModel(n_states, n_actions),
	transitions(n_states, n_actions, 0.5),
    mean_mdp(n_states, n_actions, NULL),
    reward_family(reward_family_)
{
    printf("Creating DiscreteMDPCounts with %d states and %d actions\n",  n_states, n_actions);
    N = n_states * n_actions;
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


DiscreteMDPCounts::DiscreteMDPCounts(const DiscreteMDPCounts& model) :
	MDPModel(model.n_states, model.n_actions),
	use_sampling(model.use_sampling),
	transitions(model.transitions),
	mean_mdp(model.mean_mdp),
	reward_family(model.reward_family),
	N(model.N)
{
	printf("Copying DiscreteMDPCounts with %d states and %d actions\n",  n_states, n_actions);
	
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

/// CHECK: Some parameters are not copied
DiscreteMDPCounts* DiscreteMDPCounts::Clone ()
{
	DiscreteMDPCounts* clone = new DiscreteMDPCounts(*this);
	return clone;
}



DiscreteMDPCounts::~DiscreteMDPCounts()
{
    for (int i=0; i<N; ++i) {
        delete ER[i];
    }
    //printf ("COUNTS MODEL\n");
    //ShowModel();
}

#if 0
/// Copy the mean MDP
DiscreteMDP* DiscreteMDPCounts::CreateMDP() const
{
    mdp_dbg("Making a DiscreteMDP with %d states, %d actions from model\n", n_states, n_actions);
	DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);
	CopyMeanMDP(mdp);
    return mdp;
}
#endif

void DiscreteMDPCounts::setFixedRewards(const Matrix& rewards)
{
	//logmsg("Setting fixed rewards\n");
    for (int s=0; s<n_states; ++s) {
        for (int a=0; a<n_actions; ++a)  {
            int ID = getID(s, a);
            delete ER[ID];
            ER[ID] = new UnknownSingularDistribution();
            ER[ID]->Observe(rewards(s,a));
            mean_mdp.reward_distribution.setFixedReward(s, a, rewards(s,a));
			//printf("R: %d %d %f -> %f\n",
			//	   s, a, rewards(s,a), ER[ID]->getMean());
        }
    }
}

void DiscreteMDPCounts::AddTransition(int s, int a, real r, int s2)
{
    int ID = getID (s, a);
    //printf ("(%d, %d) [%.2f] -> %d\n", s, a, r, s2);
    transitions.Observe(s, a, s2);
    ER[ID]->Observe(r);

    real expected_reward = getExpectedReward(s,a);
    mean_mdp.reward_distribution.setFixedReward(s, a, expected_reward);
    for (int s_next=0; s_next<n_states; s_next++) {
		real p = transitions.marginal_pdf(s, a, s_next);
        mean_mdp.setTransitionProbability(s, a, s_next, p);
    }
    
}

//void DiscreteMDPCounts::SetNextReward(int s, int a, real r)
//{
//    ER[getID (s, a)].mean = r;
//}

/// Generate a reward from the marginal distribution 
real DiscreteMDPCounts::GenerateReward (int s, int a) const
{
    return ER[getID (s, a)]->generateMarginal();
}

/// Generate a transition from the marginal distribution
int DiscreteMDPCounts::GenerateTransition (int s, int a) const
{
	return transitions.marginal_generate(s, a);
}

/// Get the specific transition probability
real DiscreteMDPCounts::getTransitionProbability (int s, int a, int s2) const
{
    return transitions.marginal_pdf(s, a, s2);
}

/// Get a vector of transition probabilities
Vector DiscreteMDPCounts::getTransitionProbabilities (int s, int a) const
{
    return transitions.getMarginal(s, a);
}

/// get the expected reward
real DiscreteMDPCounts::getExpectedReward (int s, int a) const
{
    return ER[getID (s,a)]->getMean();
}

/// Reset at the end of an episode 
void DiscreteMDPCounts::Reset()
{
}

/// Show the model.
void DiscreteMDPCounts::ShowModel() const
{
	printf ("# mean model\n");
    for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << "P: " << a << "," << i << ":";
            for (int j=0; j<n_states; j++) {
                real p = getTransitionProbability(i, a, j);
                //if (p<0.01) p =0.0f;
                std::cout << p << " ";
            }
            std::cout << " ["
                      << transitions.getParameters(i, a).Sum()
                      << "]\n";
        }
    }

	for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << "R(" << a << "," << i 
                      << ") = " << getExpectedReward(i, a)
				//<< " [" << ER[getID(i,a)].n_samples << "]"
                      << std::endl; 
        }
	}
}

/// Show the model.
void DiscreteMDPCounts::ShowModelStatistics() const
{
	printf ("# model statistics\n");
    for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << "P: " << a << "," << i << ":";
			transitions.getParameters(i, a).print(stdout);
        }
    }

	for (int a=0; a<n_actions; a++) {
        for (int i=0; i<n_states; i++) {
            std::cout << "R(" << a << "," << i 
                      << ") = " << getExpectedReward(i, a)
				//<< " [" << ER[getID(i,a)].n_samples << "]"
                      << std::endl; 
        }
	}
}


DiscreteMDP* DiscreteMDPCounts::generate() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions, NULL);
    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            //Vector C =  P[getID (s,a)].getMarginal();
            Vector C =  transitions.generate(s,a);
            real expected_reward = GenerateReward(s,a);
            mdp->reward_distribution.addFixedReward(s, a, expected_reward);
            for (int s2=0; s2<n_states; s2++) {
                if (C[s2]) {
                    mdp->setTransitionProbability(s, a, s2, C[s2]);
                }
            }
        }
    }
    
    return mdp;
}


/// Get a pointer to the mean MDP
const DiscreteMDP * const DiscreteMDPCounts::getMeanMDP() const
{
	//DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);
	//CopyMeanMDP(mdp);
    //    return mdp;
    return &mean_mdp;
}

void DiscreteMDPCounts::CopyMeanMDP(DiscreteMDP* mdp) const
{
    if (mdp->getNStates() != n_states) {
        throw std::runtime_error("incorrect number of states");
    }

    if (mdp->getNActions() != n_actions) {
        throw std::runtime_error("incorrect number of actions");
    }

    for (int s=0; s<n_states; s++) {
        for (int a=0; a<n_actions; a++) {
            Vector C =  transitions.getMarginal(s, a);
            real expected_reward = getExpectedReward(s,a);
            mdp->reward_distribution.addFixedReward(s, a, expected_reward);
            for (int s2=0; s2<n_states; s2++) {
                mdp->setTransitionProbability(s, a, s2, C[s2]);
            }
        }
    }
    
}


