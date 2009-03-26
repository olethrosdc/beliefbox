// -*- Mode: C++; -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BANDIT_UCT_PROBLEM_H
#define BANDIT_UCT_PROBLEM_H

#include "PolicyEvaluation.h"
#include "ValueIteration.h"
#include "BetaDistribution.h"
#include "NormalDistribution.h"
#include "Random.h"
#include "EasyClock.h"

#include <list>
#include <vector>
#include <set>

/// A prior 
class ContextBanditPrior {
protected:
    real tau;
public:
    std::vector<NormalDistributionUnknownMean> prior;
    ContextBanditPrior()
    {
    }
    ContextBanditPrior(std::vector<NormalDistributionUnknownMean> prior_)
        : prior(prior_)
    {}

    ContextBanditPrior(int n_states, int n_actions, real tau_, real mu_0, real tau_0)
    {
        tau = tau_;
        for (int i=0; i<n_states; i++) {
            for (int j=0; j<n_actions; j++) {
                //prior.push_back(BetaDistribution(alpha, beta));
                prior.push_back(NormalDistributionUnknownMean(mu_0, tau_0, tau));
            }
        }
    }
    
    NormalDistributionUnknownMean getStateActionReward(int state, int action)
    {
        assert(action >= 0 && action < (int) prior.size());
        return prior[action];
    }
        NormalDistributionUnknownMean& StateActionReward(int state, int action)
    {
        assert(action >= 0 && action < (int) prior.size());
        return prior[action];
    }
};

/// a prior distribution on the number of arms
class ContextBanditBelief
{
protected:
    int n_states, n_actions;
    ContextBanditPrior prior;
    real _greedy_return;
    bool _greedy_return_is_set;
public:
    ContextBanditBelief()
    {
        n_actions = 0;
        _greedy_return_is_set = false;
    }
    /// Create a belief
    ContextBanditBelief(int n_states_, int n_actions_, real tau, real mu_0, real tau_0)
        : n_states(n_states_),
          n_actions(n_actions_),
          prior(n_states_, n_actions_, tau, mu_0, tau_0),
          _greedy_return_is_set(false)
    {
    }
    ~ContextBanditBelief()
    {
        //for (int i=0; i<n_actions; i++) {
            //    delete prior[i];
        //}
    }
    ContextBanditPrior getPrior()
    {
        assert(n_states > 0 && n_actions > 0);
        return prior;
    }
    void update(int state, int action, real reward, int next_state)
    {
        assert (state >= 0 && state < n_states);
        assert (action >= 0 && action < n_actions);
        //assert (approx_eq(reward, 0.0) || approx_eq(reward, 1.0));
        prior.StateActionReward(state,action).calculatePosterior(reward);
        _greedy_return_is_set = false;
    }

    real getProbability(int state, int action, real reward, int next_state)
    {
        assert (action >= 0 && action < n_actions);
        real p = prior.StateActionReward(state, action).getMean();
        if (approx_eq(reward, 0.0)) {
            return 1.0 - p;
        } else if (approx_eq(reward, 1.0)) {
            return p;
        } else {
            fprintf (stderr, "Reward of %f is illegal for s,a=0,0\n", reward);
            exit(-1);
        }
        return 0.0;
    }
    

    /// Gives the reward starting from a particular state 
    /// and playing action a
    real getMeanReward(int state, int action)
    {
        return prior.StateActionReward(state, action).getMean();
    }

    /// Gives the greedy return starting from a particular statea
    real getGreedyReturn(int state, real gamma)
    {
        if (_greedy_return_is_set) {
            return _greedy_return;
        }
        std::vector<real> p(n_actions);
        int arg_max = 0;
        for (int i=0; i<n_actions; i++) {
            p[i] = prior.StateActionReward(state, i).getMean();
            if (p[i] > p[arg_max]) {
                arg_max = i;
            }
        }
        real p_max = p[arg_max];
        real R = p_max;
        _greedy_return = R / (1.0 - gamma);
        _greedy_return_is_set = true;
        return _greedy_return;
    }

    /// Gives the Thompson return
    real sampleReturn(int state, real gamma)
    {
        std::vector<real> p(n_actions);
        int arg_max = 0;
        for (int i=0; i<n_actions; i++) {
            p[i] = prior.StateActionReward(state, i).generate();
            if (p[i] > p[arg_max]) {
                arg_max = i;
            }
        }
        real p_max = p[arg_max];
        real R = p_max;
        real U = R / (1.0 - gamma);
        return U;
    }
};

template <typename ContextBanditBelief>
class BeliefTree
{
protected:
    std::vector<Distribution*> densities; ///< hold densities for the tree-derived MDPs.
public:
    class Edge;

    /// Node class
    ///
    /// Contains a list of edges and an incoming edge.
    /// It summarises the total reward and probability.
    /// It also contains upper and lower bounds on the future return.
    class Node
    {
    public:
        ContextBanditBelief belief;
        int state;
        std::list<Edge*> outs;
        Edge* in_edge;
        int index;
        int depth;
        std::vector<real> U; ///< upper bound samples
        real R; ///< reward received along path
        real U_c; ///< current upper bound
        real L; ///< lower bound
        real p; ///< probability of reaching node
        /// Constructor only sets up the probability
        Node() : p(1.0)
        {
        }
        Node* GetParent() {
            if (in_edge) {
                return in_edge->src;
            }
            return NULL;
        }
        real GetPathProbability()
        {
            return p;
        }
        real GetIncomingReward()
        {
            if (in_edge) {
                return in_edge->r;
            } else {
                return 0.0;
            }
        }
    };
    
    /// Edge class
    /// 
    /// Information about links between nodes.
    /// Also stores the observations used to get from
    /// the previous node to the current.
    class Edge
    {
    public:
        Node* src; ///< source node
        Node* dst; ///< destination node
        int a; ///< action taken
        real r; ///< reward received
        real p; ///< probability of path component
        real GetEdgeProbability()
        {
            return p;
        }
        Edge (Node* src_, Node* dst_, int a_, real r_, real p_)
            : src(src_), dst(dst_), a(a_), r(r_), p(p_)
        {
        }
    };

    Node* root;
    
    std::list<Node*> nodes; ///< a list of nodes for book-keeping purposes
    //std::list<Edge*> edges; ///< a list of edges for book-keeping purposes
    
    int n_states;
    int n_actions;
    real gamma;
    BeliefTree(ContextBanditBelief prior,
               int state,
               int n_states_,
               int n_actions_,
               real gamma_)
        :
        n_states(n_states_),
        n_actions(n_actions_),
        gamma(gamma_)
    {
        root = new Node;
        root->belief = prior;
        root->state = state;
        root->index = 0;
        root->depth = 0;
        root->in_edge = NULL;

        nodes.push_back(root);
    }

    ~BeliefTree()
    {
        DeleteDensities();

        for (typename std::list<Node*>::iterator i=nodes.begin();
             i!=nodes.end(); ++i) {
            Node* n = *i;
            DeleteNode(n);
        }
    }

    // Should be called after created MDPs are discarded
    void DeleteDensities()
    {
        for (std::vector<Distribution*>::iterator i=densities.begin();
             i!=densities.end();
             ++i) {
            Distribution *d = *i;
            delete d;
        }
    }

    int GetNumberOfLeafNodes()
    {
        int n_leaf_nodes = 0;
        // leaf ndoes are state nodes
        for (typename std::list<Node*>::iterator i=nodes.begin();
             i!=nodes.end(); ++i) {
            Node* node = *i;
            if (node->outs.size() == 0) {
                n_leaf_nodes++;
            }
        }
        return n_leaf_nodes;
    }

    /// Return the newly created node
    Node* ExpandAction(Node *selected_node, int a, real r, int s, int verbose = 0)
    {
        Node* next = new Node;
        next->belief = selected_node->belief;
        next->state = s;
        next->index = nodes.size();
        next->depth = selected_node->depth + 1;
        next->R = selected_node->R + r * pow(gamma, next->depth);
        // the probability of the next state and reward given the
        // belief, state and action
        real p = selected_node->belief.getProbability(selected_node->state, a, r, s);
        real p_path = p * selected_node->GetPathProbability();
        next->p = p_path; // fill in the total path probability
        next->belief.update(selected_node->state, a, r, s); // update the belif
        
        // save the edge connecting the previous node to the next
        Edge* next_edge = new Edge(selected_node,
                                 next, //was: nodes[nodes.size()-1],
                                 a, // action taken
                                 r, // reward observed
                                 p // probability given previous node and action
                                   );
        //edges.push_back(next_edge);


        if (verbose >= 100) {
            printf ("Added edge  %d --(%d %f %f)-> %d\n",
                    next_edge->src->index,
                    next_edge->a,
                    next_edge->r,
                    next_edge->p,
                    next_edge->dst->index);
        }
                
        // save the edge to the list of output edges of the previous node
        // and as an input edge of the next node
        selected_node->outs.push_back(next_edge);
        next->in_edge = next_edge;
        
        // save the node
        nodes.push_back(next);

        // return it
        return nodes.back();
    }

    /// Find the next node
    Node* FindObservation(Node* src, int a, real r, int s, int verbose = 0)
    {
        for (typename std::list<Edge*>::iterator j = src->outs.begin();
             j != src->outs.end(); ++j) {
            Edge* edge = *j;
            if (edge->a == a
                && edge->r == r
                && edge->dst->state == s) {
                return edge->dst;
            }
        }
        
        return NULL;
    }
        /// Cut a tree, making node i the root
    void MakeRoot(Node* node, int verbose = 0)
    {
        Node* parent = node->GetParent();
        if (!parent) {
            if (verbose >= 50) {
                printf("# warning: making root a node which is already root\n");
            }
            return; // node is already the root node
        }
        

        RecursiveDeleteExcept(root, node);
        root = node;
        root->index = 0;
        root->depth = 0;
        root->in_edge = NULL;
        UpdateNodes(root, 0);
    }
    
    /// Update index and depth information of nodes
    ///
    /// Fills out index numbers in a depth-first manner.
    /// This if index1 > index2, depth1 >= depth2
    int UpdateNodes(Node* node, int index)
    {
        int depth = node->depth;
        for (typename std::list<Edge*>::iterator j = node->outs.begin();
             j != node->outs.end(); ++j) {
            Edge* edge = *j;
            Node* next = edge->dst;
            next->index = ++index;
            next->depth = depth + 1;
        }

        for (typename std::list<Edge*>::iterator j = node->outs.begin();
             j != node->outs.end(); ++j) {
            Edge* edge = *j;
            Node* next = edge->dst;
            index = UpdateNodes(next, index);
        }
        return index;
    }

    /// Delete a node and outgoing edges
    ///
    /// Assumes later nodes are deleted first,
    /// otherwise edge data is lost.
    void DeleteNode(Node* node)
    {
        
        assert(node);
        for (typename std::list<Edge*>::iterator j = node->outs.begin();
             j != node->outs.end(); ++j) {
            Edge* edge = *j;
            delete edge;
        }
        delete node;
    }

    void RecursiveDeleteExcept(Node* node, Node* exception)
    {
        if (node==exception) {
            return;
        }
        for (typename std::list<Edge*>::iterator j = node->outs.begin();
             j != node->outs.end(); ++j) {
            Edge* edge = *j;
            Node* next = edge->dst;
            if (next != NULL && next != exception) {
                RecursiveDeleteExcept(next, exception);
            }
        }
        DeleteNode(node);
        nodes.remove(node);
    }
    /// Expand a node in the tree
    void Expand(Node* node, int verbose = 0)
    {
        for (int a=0; a<n_actions; a++) {
            // If we play, there are two possibilities
            ExpandAction(node, a, 1.0, 0, verbose);
            ExpandAction(node, a, 0.0, 0, verbose); 
        }
    }

    int FindRootAction(Node* node)
    {
        int a = -1;
        while (node->in_edge) {
            a = node->in_edge->a;
            if (node->in_edge->src) {
                node = node->in_edge->src;
            } else {
                assert(node->in_edge->src == root);
            }
        }
        return a;
    }

    std::list<Node*>& getNodes()
    {
        return nodes;
    }

    DiscreteMDP CreateMeanMDP(real gamma, int verbose = 0)
    {
        int n_nodes = nodes.size();
        int terminal = n_nodes;
		
        // clear mean MDP
        DiscreteMDP mdp(n_nodes + 1, n_actions, NULL, NULL);

        // no reward in the first state
        {
            Distribution* reward_density = 
                new SingularDistribution(0.0);
            densities.push_back(reward_density);
            for (int a=0; a<n_actions; a++) {
                mdp.setRewardDistribution(0,
                                          a,
                                          reward_density);
            }
        }

        for (typename std::list<Node*>::iterator i=nodes.begin(); i!=nodes.end(); ++i) {
            Node* node = *i;
            int n_edges = node->outs.size();
            if (verbose >= 90) {
                printf ("Node %d has %d outgoing edges\n",
                        node->index, n_edges);
            }
            // loop for internal nodes
            for (typename std::list<Edge*>::iterator j = node->outs.begin();
                 j != node->outs.end(); ++j) {
                Edge* edge = *j;
                Distribution* reward_density = 
                    new SingularDistribution(edge->r);
                
                densities.push_back(reward_density);
                mdp.setTransitionProbability(node->index,
                                             edge->a,
                                             edge->dst->index,
                                             edge->p);
                if (verbose >= 90) {
                    printf ("%d: - a=%d - r=%f - p=%f ->%d\n", node->index, edge->a, edge->r, edge->p, edge->dst->index);
                }
                for (int a=0; a<n_actions; a++) {
                    mdp.setRewardDistribution(edge->dst->index,
                                              a,
                                              reward_density);
                }
            }
           
            // the leaf nodes
            if (!n_edges) {
                real mean_return = node->belief.getGreedyReturn(node->state, gamma);
                
                for (int a=0; a<n_actions; a++) {
                    mdp.setTransitionProbability(node->index, a, terminal, 1.0);
                    // the actions are fake, we only use the _previous_ reward!
                    real r = node->GetIncomingReward();
                    Distribution* reward_density = 
                        new SingularDistribution(r + gamma*mean_return);
                    //new SingularDistribution(mean_reward);
                    densities.push_back(reward_density);
                    mdp.setRewardDistribution(node->index, a, reward_density);
                }
            }
        }


        // the terminal node
        for (int a=0; a<n_actions; a++) {
            Distribution* reward_density = 
                new SingularDistribution(0.0);
            densities.push_back(reward_density);
            mdp.setTransitionProbability(terminal, a, terminal, 1.0);
            mdp.setRewardDistribution(terminal, a, reward_density);
        }
        
        return mdp;
    }



    DiscreteMDP CreateUpperBoundMDP(real gamma, int verbose = 0)
    {
        int n_nodes = nodes.size();
        int terminal = n_nodes;

        if (verbose >= 90) {
            printf ("Creating MDP with %d nodes\n", n_nodes);
        }
        DiscreteMDP mdp(n_nodes + 1, n_actions, NULL, NULL);
        // assume MDP is cleared
        // no reward in the first state
        {
            Distribution* reward_density = 
                new SingularDistribution(0.0);
            densities.push_back(reward_density);
            for (int a=0; a<n_actions; a++) {
                mdp.setRewardDistribution(0,
                                          a,
                                          reward_density);
            }
        }

        for (typename std::list<Node*>::iterator i = nodes.begin();
             i!=nodes.end(); ++i) {
            Node* node = *i;
            int n_edges = node->outs.size();
            if (verbose >= 90) {
                printf ("Node %d has %d outgoing edges\n", node->index, n_edges);
            }
            for (typename std::list<Edge*>::iterator j = node->outs.begin();
                 j != node->outs.end(); ++j) {
                Edge* edge = *j;
                Distribution* reward_density = 
                    new SingularDistribution(edge->r);
                
                densities.push_back(reward_density);
                mdp.setTransitionProbability(node->index,
                                             edge->a,
                                             edge->dst->index,
                                             edge->p);
                if (verbose >= 90) {
                    printf ("%d: - a=%d - r=%f - p=%f ->%d\n", node->index, edge->a, edge->r, edge->p, edge->dst->index);
                }
                for (int a=0; a<n_actions; a++) {
                    mdp.setRewardDistribution(edge->dst->index,
                                              a,
                                              reward_density);
                }
            }
            
            if (!n_edges) {
                real mean_reward = node->belief.getGreedyReturn(node->state, gamma);
                while (node->U.size() <= 0) {
                    node->U.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ub = Mean(node->U);
                if (Ub < mean_reward) {
                    Ub = mean_reward;
                }

                for (int a=0; a<n_actions; a++) {
                    mdp.setTransitionProbability(node->index, a, terminal, 1.0);
                    real r = node->GetIncomingReward();
                    Distribution* reward_density = 
                        new SingularDistribution(r + gamma * Ub);
                    densities.push_back(reward_density);
                    mdp.setRewardDistribution(node->index, a, reward_density);
                }
            }
        }

        for (int a=0; a<n_actions; a++) {
            Distribution* reward_density = 
                new SingularDistribution(0.0);
            densities.push_back(reward_density);
            mdp.setTransitionProbability(terminal, a, terminal, 1.0);
            mdp.setRewardDistribution(terminal, a, reward_density);
        }
        
        return mdp;
    }

};


enum ExpansionMethod {
    SerialExpansion = 0, // (54.9)
    RandomExpansion,  //1
    HighestMeanValue, //2 (53.5)
    HighestDiscountedMeanValue, //3
    ThompsonSampling, // 4
    DiscountedThompsonSampling, //5
    ThompsonBound, // 6
    DiscountedThompsonBound, //7
    HighProbabilityBound, // 8 (55.5)
    DiscountedHighProbabilityBound, // 9 (53.8)
    HighProbabilityBoundOnly, // 10 (54.0)
    DiscountedHighProbabilityBoundOnly, // 11
    MeanHighProbabilityBound, // 12
    DiscountedMeanHighProbabilityBound, // 13
    GreedyBoundReduction, // 14 (54.3)
    BAST, // 15 
    BAST_FULL // 16
};

int MakeDecision(BeliefTree<ContextBanditBelief>& new_tree,
                 ExpansionMethod expansion_method,
                 int n_states,
                 int n_actions,
                 ContextBanditBelief prior,
                 int state,
                 real gamma,
                 int n_iter,
                 int verbose,
                 int max_value_iterations,
                 FILE* fout = NULL);

typedef BeliefTree<ContextBanditBelief>::Node BeliefTreeNode;
typedef BeliefTree<ContextBanditBelief>::Edge BeliefTreeEdge;
typedef std::list<BeliefTreeNode*> BTNodeSet;
typedef std::list<BeliefTreeEdge*> BTEdgeSet;

#endif
