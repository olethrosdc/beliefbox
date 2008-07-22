// -*- Mode: C++; -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "PolicyEvaluation.h"
#include "BetaDistribution.h"
#include "Random.h"

#include <list>
#include <vector>
#include <set>

class SimpleBelief
{
protected:
    BetaDistribution prior;
    real r1;
    real r2;
public:
    /// Create a belief
    SimpleBelief(real alpha = 1.0, real beta=1.0, real r1_=0.0, real r2_=1.0)
        : prior(alpha, beta), r1(r1_), r2(r2_)
    {
    }

    void update(int state, int action, real reward, int next_state)
    {
        if (state==0 && action == 0) {
            if (reward == r1) {
                prior.calculatePosterior(0.0);
            } else if (reward == r2) {
                prior.calculatePosterior(1.0);
            } else {
                fprintf (stderr, "Reward of %f should not have been observed in this state\n", reward);
            }
        }
    }

    real getProbability(int state, int action, real reward, int next_state)
    {
        if (state==0 && action == 0) {
            real p = prior.getMean();
            if (reward == r1) {
                return 1.0 - p;
            } else if (reward == r2) {
                return p;
            } else {
                fprintf (stderr, "Reward of %f is illegal for s,a=0,0\n", reward);
            }
        }
        return 1.0;
    }
    
};

template <typename B>
class BeliefTree
{
public:
    class Edge;
    class Node
    {
    public:
        B belief;
        int state;
        std::list<Edge*> edges;
    };
    
    class Edge
    {
    public:
        Node* src;
        Node* dst;
        real r;
        real p;
    };

    std::vector<Node> nodes;
    std::vector<Edge> edges;
    
    int n_states;
    int n_actions;
    BeliefTree(SimpleBelief prior,
               int state,
               int n_states_,
               int n_actions_)
        :
        n_states(n_states_),
        n_actions(n_actions_)
    {
        Node root;
        root.belief = prior;
        root.state = state;
        nodes.push_back(root);
    }
    /// Return 
    Node& ExpandAction(int i, int a, real r, int s)
    {
        Node next;
        next.belief = nodes[i].belief;
        next.state = s;
        real p = nodes[i].belief.getProbability(nodes[i].s, a, r, s);
        next.belief.update(nodes[i].s, a, r, s);
        Node& node_ref = nodes.push_back(next);
        Edge edge;
        edge.src = nodes[i];
        edge.dst = node_ref;
        edge.p = p;
        Edge edge_ref=  edges.push_back(edge);
        nodes[i].edges.push_back(edge_ref);
        return node_ref;
    }

    /// Return 
    std::set<Node&> Expand(int i)
    {
        B belief = nodes[i].belief;
        int state = nodes[i].state;
        std::set<Node&> node_set;

        if (state==1) { 
            // terminal state, all actions do nothing whatsoever.            
            node_set.push_back(ExpandAction(i, 0, 0.0, 1));
            node_set.push_back(ExpandAction(i, 1, 0.0, 1));
        } else if (state==0) {
            // If we play, there are two possibilities
            node_set.push_back(ExpandAction(i, 0, 1.0, 0));
            node_set.push_back(ExpandAction(i, 0, -1.0, 0));
            // if we move to the terminal state nothing happens
            node_set.push_back(ExpandAction(i, 1, 0.0, 1));
        }
        return node_set;
    }


    DiscreteMDP CreateMDP()
    {
        int current_node = nodes[0];
        int n_nodes = nodes.size();
        DiscreteMDP (n_nodes, 2, NULL);
    }
};



/// A toy UCT stopping problem




//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int main (int argc, char** argv)
{
    real alpha = 0;
    real beta = 0;
    real actual_probability = urandom(0, 1);
    
    SimpleBelief prior(1, 1, -1.0, 1.0);

    BeliefTree<SimpleBelief> test(prior, 0, 2, 2);
    return 0;
}


#endif
