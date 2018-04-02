// -*- Mode: c++ -*-
// copyright (c) 2017 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "TreeBRL.h"


//#define TBRL_DEBUG

TreeBRL::TreeBRL(int n_states_,
                 int n_actions_,
                 real gamma_,
                 MDPModel* belief_,
                 RandomNumberGenerator* rng_,
                 int horizon_,
				 enum LeafNodeValue leaf_node)
    : n_states(n_states_),
      n_actions(n_actions_),
      gamma(gamma_),
      belief(belief_),
      rng(rng_),
      horizon(horizon_),
      T(0),
      size(0),
      Qs(n_actions),
	  leaf_node_expansion(leaf_node)
{
	const char* leaf_value_name[] = {"None", "V_Min", "V_Max", "V_mean", "V_U", "V_L"};
    logmsg("Starting Tree-Bayes-RL with %d states, %d actions, %d horizon, %s bounds\n", n_states, n_actions, horizon, leaf_value_name[leaf_node]);

    current_state = -1;

}

// Note that the belief tree is only created within Act() and
// destroyed immediately. Hence there is no need to remove
// anything here. But on the other hand, this is inefficient.
TreeBRL::~TreeBRL()
{
    //printf(" # destroying tree of size %d\n", size);
}

void TreeBRL::Reset()
{
    current_state = -1;
    current_action = -1;
}

void TreeBRL::Reset(int state)
{
    current_state = state;
    current_action = -1;
}

/// Full observation
real TreeBRL::Observe (int state, int action, real reward, int next_state, int next_action)
{
    if (state>=0) {
        belief->AddTransition(state, action, reward, next_state);
    }
    current_state = next_state;
    current_action = next_action;
    return 0.0;
}
/// Partial observation 
real TreeBRL::Observe (real reward, int next_state, int next_action)
{
    if (current_state >= 0 && current_action >= 0) {
        belief->AddTransition(current_state, current_action, reward, next_state);
    }
    current_state = next_state;
    current_action = next_action;
    return 0.0;
}


/// Get an action using the current exploration policy.
/// it calls Observe as a side-effect.
int TreeBRL::Act(real reward, int next_state)
{
    assert(next_state >= 0 && next_state < n_states);

    T++;
    if (current_state >= 0 && current_action >= 0) {
        belief->AddTransition(current_state, current_action, reward, next_state);
    }

    current_state = next_state;
    
    //int n_MDP_leaf_samples = 1;
    BeliefState belief_state = CalculateBeliefTree();
	
	//printf("%f %f %f\n", belief_state.CalculateValues(), 
	//belief_state.CalculateValues(leaf_node_expansion);

	//Qs.printf(stdout);
    int next_action = ArgMax(Qs);
	//printf("-> %d\n", next_action);
	// sometimes act randomly
	if (rng->uniform() < 0) {
		next_action = rng->random() % n_actions;
	}
    current_action = next_action;
    return current_action;
}


    /// Calculate a sparse belief tree where we take n_samples state
    /// samples and use n_TS MDP samples for the upper and lower bounds at
    /// the leaf nodes

void TreeBRL::CalculateSparseBeliefTree(int n_samples, int n_TS)
{
    // Initialise the root belief state
    BeliefState belief_state(*this, belief, current_state);
    belief_state.SparseExpandAllActions(n_samples);
    belief_state.CalculateValues(leaf_node_expansion);
    //belief_state.CalculateLowerBoundValues(n_TS),
    //belief_state.CalculateUpperBoundValues(n_TS));
}

/// Calculate a belief tree.
TreeBRL::BeliefState TreeBRL::CalculateBeliefTree()
{
    // Initialise the root belief state
    BeliefState belief_state(*this, belief, current_state);
    belief_state.ExpandAllActions();
	belief_state.CalculateValues(leaf_node_expansion);
	return belief_state;
}

//------------- Belief states ----------------//

TreeBRL::BeliefState::BeliefState(TreeBRL& tree_,
                                  MDPModel* belief_,
                                  int state_) : tree(tree_), state(state_), probability(1), t(0)
{
	belief = belief_->Clone(); 
    tree.size++;
}

/// Use this to construct a subsequent belief state
TreeBRL::BeliefState::BeliefState(TreeBRL& tree_,
                                  MDPModel* belief_,
                                  int prev_state_,
                                  int prev_action_,
                                  int state_,
                                  real r,
                                  real p,
                                  BeliefState* prev_)
	: tree(tree_), 
	  state(state_),
	  prev_action(prev_action_),
	  prev_reward(r), probability(p), prev(prev_), t(prev_->t + 1)
{
	
#ifdef TBRL_DEBUG
	logmsg("Cloning belief");
#endif
	belief = belief_->Clone();
#ifdef TBRL_DEBUG
	logmsg("Adding new transition");
#endif
    belief->AddTransition(prev_state_,
						  prev_action,
						  prev_reward,
						  state);
	tree.size++;
#ifdef TBRL_DEBUG
	logmsg("Previous belief");
	belief_->ShowModelStatistics();
	logmsg("Next belief");
	belief->ShowModelStatistics();
	logmsg("Tree size: %d\n", tree.size);
#endif

}

TreeBRL::BeliefState::~BeliefState()
{
	tree.size--;
	//delete belief;
}


/// Generate transitions from the current state for all
/// actions. Do this recursively, using the marginal
/// distribution, but using sparse sampling.
///
void TreeBRL::BeliefState::SparseExpandAllActions(int n_samples)
{
    if (t >= tree.horizon) {
        return;
    }
    real p = 1 / (real) n_samples;
    for (int k=0; k<n_samples; ++k) {
        for (int a=0; a<tree.n_actions; ++a) {
            int next_state = belief->GenerateTransition(state, a);
            real reward = belief->GenerateReward(state, a);
            // Generate the new belief state and put it in the tree
            children.push_back(BeliefState(tree, belief, state, a, next_state, reward, p, this));
        }
    }

    for (uint i=0; i<children.size(); ++i) {
        children[i].SparseExpandAllActions(n_samples);
    }
}
/// Generate transitions from the current state for all
/// actions. Do this recursively, using the marginal distribution. Only expand the children when we're under the horizon.
void TreeBRL::BeliefState::ExpandAllActions()
{
    if (t >= tree.horizon) {
        return;
    }
    for (int a=0; a<tree.n_actions; ++a) {
        for (int next_state=0;
             next_state<tree.n_states;
             ++next_state) {
            real p = belief->getTransitionProbability(state, a, next_state);
            real reward = belief->GenerateReward(state, a);
            children.push_back(BeliefState(tree, belief, state, a, next_state, reward, p, this));
        }
    }
            
    for (uint i=0; i<children.size(); ++i) {
        children[i].ExpandAllActions();
    }
}

            
        
/// Return the values using Backwards induction on the already
/// constructed MDP.
///
/// If w is the current belief state, and w' the successor, while a is
/// our action, then we can write the value function recursion:
///
/// \f$V_t(w) = \max_a Q_t(w, a) = \max_a E\{r(w,a,w') + \gamma V_{t+1} (w')\}\f$,
///
/// where the expectation is 
/// \f$Q_t(w, a) = \sum_{s'} {r(w,a,s') + \gamma P(s' | a, s) V_{t+1} (w')\}\f$ and \f$w' = w( | s, a, s')\f$.
real TreeBRL::BeliefState::CalculateValues(LeafNodeValue leaf_node)
{
    Vector Q(tree.n_actions);
    real V = 0;
    real discount = tree.gamma;
	
    if (t < tree.horizon) {
        for (uint i=0; i<children.size(); ++i) {
            int a = children[i].prev_action;
			real p = children[i].probability;
			real r = children[i].prev_reward;
			int s_next = children[i].state;
			real V_next = children[i].CalculateValues(leaf_node);
            Q(a) += p * (r + discount * V_next);
#ifdef TBRL_DEBUG
			printf("t:%d s:%d i:%d a:%d p:%f s2:%d, r:%f v:%f\n",
				   t, state, i, a, p, s_next, r, V_next);
#endif
        }
        V += Max(Q);
#ifdef TBRL_DEBUG
		Q.print(stdout); printf(" %d/%d\n", t, tree.horizon);
#endif
    } else {
		switch(leaf_node) {
		case NONE: V = 0; break;
		case V_MIN: V = 0; break;
		case V_MAX: V = 1.0 / (1.0 - discount); break;
		case V_MEAN: V = MeanMDPValue(); break;
		case V_UTS: V = UTSValue(); break;
		case V_LTS: V = LTSValue(); break;
		}
    }
	


    if (t==0) {
        tree.Qs = Q;
    }
    return V;
}

/// Return the values using an upper bound
real TreeBRL::BeliefState::UTSValue()
{
	int n_samples = 2;
	real V = 0;
	for (int i=0; i<n_samples; ++i) {
		DiscreteMDP* model = belief->generate();
		ValueIteration VI(model, tree.gamma);
		VI.ComputeStateValuesStandard(1e-3);
		V += VI.getValue(state);
		delete model;
	}
	V /=  (real) n_samples;
	return V;
}

/// Return the values using a lower bound
real TreeBRL::BeliefState::LTSValue()
{
    real discount = tree.gamma;
	int n_samples = 2;
	const DiscreteMDP* model = belief->getMeanMDP();
	ValueIteration VI(model, discount);
	VI.ComputeStateValuesStandard(1e-3);
	FixedDiscretePolicy* policy = VI.getPolicy();
	real V_next = 0;
	for (int i=0; i<n_samples; ++i) {
		DiscreteMDP* model = belief->generate();
		PolicyEvaluation PI(policy, model, tree.gamma);
		PI.ComputeStateValues(1e-3);
		V_next += PI.getValue(state);
		delete model;
	}
	real V = V_next / (real) n_samples;

    return V;
}
        
/// Return the values using the mean MDP
real TreeBRL::BeliefState::MeanMDPValue()
{
    real discount = tree.gamma;
	const DiscreteMDP* model = belief->getMeanMDP();
	ValueIteration VI(model, discount);
	VI.ComputeStateValuesStandard(1e-3);
	return VI.getValue(state);
}
        




