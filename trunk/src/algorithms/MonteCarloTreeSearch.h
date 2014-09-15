/* -*- Mode: C++; -*- */
// copyright (c) 2014 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MONTE_CARLO_TREE_SEARCH_H
#define MONTE_CARLO_TREE_SEARCH_H

#include "MersenneTwister.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFile.h"

#include "RandomPolicy.h"

#include "real.h"
#include "Matrix.h"
#include "Grid.h"
#include "Environment.h"
#include "Random.h"
#include <limits>

template <class S, class A>
class MonteCarloTreeSearch
{
private:
  real gamma;                    // Discount factor.
  ContinuousStateEnvironment* environment; // The environment model.
  RandomNumberGenerator* rng;    // Random Number generator.
  AbstractPolicy<S, A>& policy;
  int MaxDepth;                  // Maximum tree depth.
  int NRollouts;                 // Number of sampled rollouts.
  int nActions;
public:
  struct Node {
    int depth;     // Depth
    S state;
    double reward; //Reward received during the transition from the father's state  
    Node* father;  //Contains the previous state at the trajectory.
    MonteCarloTreeSearch& tree;
    bool terminal; //Represents if the node is a terminal state.
    bool leaf;     //Represents if the specific node is a leaf or not

    double epsilon;
    
    std::vector<Node*> children; //Pointer to the childrens
    int nVisits;  // Number of visits
    double aveValue; // Mean Value

    // Constructor
    Node(const int& depth_, const S& state_, const double& reward_, MonteCarloTreeSearch::Node* const father_ , MonteCarloTreeSearch& tree_, const bool& terminal_ = false)
      :	depth(depth_),   
        state(state_),   //Node's state representation.
	reward(reward_), //Reward received during the transition from the father node.
	father(father_), //Father node.
	tree(tree_),
	terminal(terminal_)
    {
      epsilon  = 1E-06;
      aveValue = 0; 
      nVisits  = 0; 
      leaf     = true;
      children.resize(tree.nActions,NULL);
    }
    
    Node(const S& state_, MonteCarloTreeSearch& tree_) 
      : state(state_),   //Node's state representation.
	tree(tree_)
    {
      depth    = 0;
      reward   = 0;
      father   = NULL;
      terminal = false;
      aveValue = 0;
      nVisits  = 0;
      leaf     = true;
      children.resize(tree.nActions,NULL);
    }

    // Destructor
    ~Node() {
      for(int i=0; i < tree.nActions; ++i)
	delete children[i];
    }

    void selectAction() {
      Node* cur = this;
      double RollingValue = 0;

      // Selection phase (we cross the tree according to the UCT)
      // printf("UCT Search\n");
      while(cur->isLeaf() == false) {
	cur = cur->UCTsearch(); //Tree policy
      }

      if(cur->isTerminal() == true) {
      	RollingValue = cur->getReward();
      } else {
	// Expansion phase
	//	while(cur->getDepth()+1 < tree.MaxDepth && cur->isTerminal() == false) {
	  cur = cur->expand();
	  //} 
	RollingValue = rollOut(cur->state);
      }

      // Backpropagation phase
      while(cur != NULL) {
	cur->updateStats(RollingValue);
	RollingValue = cur->getReward() + tree.gamma*RollingValue;
	cur = cur->father;
      }	
    }

    //Tree expansion
    Node* expand() {
      std::vector<int> Unvisited; //Pointer to the childrens
      int action;
      // The unvisited children are declared
      for(action = 0; action < tree.nActions; ++action) {
	if(children[action] == NULL) {
	  Unvisited.push_back(action);
	}
      }
      //  printf("Action selection\n");
      // We select one child (uniform) randomly among the unvisited children.
      if(Unvisited.empty()) {
	real bestValue = -1000000000;
	for(int i = 0; i < tree.nActions; i++) {
	  real curValue = (children[i]->getReward() + tree.gamma*children[i]->aveValue) + 1000*sqrt((sqrt(2)*log(nVisits)) / (children[i]->nVisits)); //UCT Search
	  if(curValue > bestValue) {
	    action = i;
	    bestValue = curValue;
	  }
	}
      } else {	
	action = Unvisited[tree.rng->discrete_uniform(Unvisited.size())];
      }
      tree.environment->Reset();
      tree.environment->setState(state);
      bool running    = tree.environment->Act(action);
      S child_state   = tree.environment->getState();
      real reward     = tree.environment->getReward();

      children[action] = new Node(depth + 1, child_state, reward, this, tree, !running); // The specific child is created

      return children[action];
    }
    
    //Tree policy
    Node* UCTsearch() {
      real bestValue = -1000000000000;
      Node* selected = NULL;
      for(int i = 0; i < tree.nActions; ++i) {
	real curValue = (children[i]->getReward() + tree.gamma*children[i]->aveValue) + 1000*sqrt((log(nVisits)) / (children[i]->nVisits)); //UCT Search
	if(curValue > bestValue) {
	  selected = children[i];
	  bestValue = curValue;
	}
      }
      return selected;
    }
    double rollOut(S state_) {
      int t = 0;
      int horizon = 1000; // tree.MaxDepth - depth;
      tree.environment->Reset();
      tree.environment->setState(state_);
      tree.policy.Reset();
      bool running = true;
      real discount = 1.0;
      real total_reward = 0.0;
      real discounted_reward = 0.0;
    
      real reward = 0.0;
      do {
	//get current state
	S state = tree.environment->getState();
              
	//choose an action using Random policy
	tree.policy.Observe(reward, state);
	A action = tree.policy.SelectAction();

	//execute the selected action
	running = tree.environment->Act(action);
              
	// get reward
	reward = tree.environment->getReward();
              
	total_reward += reward;
	discounted_reward += discount * reward;
              
	discount *= tree.gamma;

	++t;              
	if (t >= horizon) {
	  running = false;
	}
      }while(running);
      //printf("Gamma = %f, discounted_reward = %f\n",tree.gamma,discounted_reward);
      return discounted_reward;
    }
    bool isTerminal() {
      return terminal;
    }
    bool isRoot() {
      if(father == NULL)
	return true;
      else 
	return false;
    }
    bool isLeaf() {
      for(int action = 0; action < tree.nActions; ++action) {
	if(children[action] == NULL) {
	  return true;
	}
      }
      return false;
    }    
    int getDepth() {
      return depth;
    }
    double getReward() {
      return reward;
    }
    double getValue() {
      return aveValue;
    }
    void setFather(Node* father_){
      father = father_;
    }
    void setState(S state_) {
      state = state_;
    }
    void updateStats(double value) {
      nVisits++;
      aveValue += (value - aveValue)/ nVisits; // Mean value
    }
  };

  //Constructor
  MonteCarloTreeSearch(const real& gamma_, ContinuousStateEnvironment* environment_, RandomNumberGenerator* rng_, AbstractPolicy<S, A>& policy_, const int& MaxDepth_ =  100, const int& NRollouts_ = 1000)
    :gamma(gamma_),
     environment(environment_),
     rng(rng_),
     policy(policy_),
     MaxDepth(MaxDepth_),
     NRollouts(NRollouts_)
  {
    nActions = environment->getNActions(); 
  };

  //Destructor
  ~MonteCarloTreeSearch(){
    delete root;
  };

  int SelectAction(S state_) {
 
    root = new Node(0, state_, 0.0, NULL, *this);
  
    for(int i=0; i<NRollouts; ++i) {
      root->selectAction();
    }

    int sel_action = 0;
    double bestValue = root->children[0]->reward + gamma*root->children[0]->aveValue;
    
    //Find the best among the available actions
    for(int action = 1; action < nActions; ++action) {
      double curValue = root->children[action]->reward + gamma*root->children[action]->aveValue;
      if(curValue > bestValue) {
	sel_action = action;
	bestValue = curValue;
      }
    }
    environment->Reset();
    environment->setState(state_);

    delete root;

    return sel_action;
  };
protected:
  std::vector< std::vector<Node*> > levels;
  Node* root;
};
#endif
