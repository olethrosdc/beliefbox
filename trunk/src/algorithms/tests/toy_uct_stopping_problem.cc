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
#include "ValueIteration.h"
#include "BetaDistribution.h"
#include "Random.h"
#include "EasyClock.h"

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
    
    real getGreedyReturn(int state, real gamma)
    {
        if (state == 0) {
            real p = prior.getMean();
            real R = p*r2 + (1.0 - p)*r1;
            real U = R / (1.0 - gamma);
            if (U > 0) {
                return U;
            } 
        } 
        return 0.0;
    }
};

template <typename B>
class BeliefTree
{
protected:
    std::vector<Distribution*> densities;
public:
    class Edge;
    class Node
    {
    public:
        B belief;
        int state;
        std::vector<Edge*> outs;
        int index;
    };
    
    class Edge
    {
    public:
        Node* src; ///< source node
        Node* dst; ///< destination node
        int a; ///< action taken
        real r; ///< reward received
        real p; ///< probability of path
        Edge (Node* src_, Node* dst_, int a_, real r_, real p_)
            : src(src_), dst(dst_), a(a_), r(r_), p(p_)
        {
        }
    };

    std::vector<Node*> nodes;
    std::vector<Edge*> edges;
    
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
        Node* root = new Node;
        root->belief = prior;
        root->state = state;
        root->index = 0;
        nodes.push_back(root);
    }

    ~BeliefTree()
    {
        //std::cout << "D: " << densities.size() << std::endl;
        for (int i=densities.size() - 1; i>=0; --i) {
            //std::cout << i << std::endl;
            delete densities[i];
        }

        for (int i=nodes.size() - 1; i>=0; --i) {
            delete nodes[i];
        }

        for (int i=edges.size() - 1; i>=0; --i) {
            delete edges[i];
        }
    }
    /// Return 
    Node* ExpandAction(int i, int a, real r, int s, int verbose = 0)
    {
        Node* next = new Node;
        next->belief = nodes[i]->belief;
        next->state = s;
        next->index = nodes.size();
        real p = nodes[i]->belief.getProbability(nodes[i]->state, a, r, s);
        next->belief.update(nodes[i]->state, a, r, s);

        
        edges.push_back(new Edge(nodes[i],
                                 next, //nodes[nodes.size()-1],
                                 a, r, p));

        int k = edges.size() - 1;
        if (verbose >= 100) {
            printf ("Added edge %d : %d --(%d %f %f)-> %d\n",
                    k,
                    edges[k]->src->index,
                    edges[k]->a,
                    edges[k]->r,
                    edges[k]->p,
                    edges[k]->dst->index);
        }
                
        nodes[i]->outs.push_back(edges[k]);

        nodes.push_back(next);
        return nodes.back();
    }

    /// Return 
    //    std::vector<Node&> Expand(int i)
    void Expand(int i, int verbose = 0)
    {
        //B belief = nodes[i]->belief;
        int state = nodes[i]->state;
        //std::vector<Node*> node_set;
        if (state==1) { 
            // terminal state, all actions do nothing whatsoever.            
            ExpandAction(i, 0, 0.0, 1, verbose);
            ExpandAction(i, 1, 0.0, 1, verbose);
        } else if (state==0) {
            // If we play, there are two possibilities
            ExpandAction(i, 0, 1.0, 0, verbose);
            ExpandAction(i, 0, -1.0, 0, verbose); // VALGRIND
            // if we move to the terminal state nothing happens
            ExpandAction(i, 1, 0.0, 1, verbose);
        }
        //return node_set;
    }

    std::vector<Node*>& getNodes()
    {
        return nodes;
    }

    DiscreteMDP CreateMDP(real gamma, int verbose = 0)
    {
        int n_nodes = nodes.size();
        int terminal = n_nodes;

        DiscreteMDP mdp(n_nodes + 1, 2, NULL, NULL);
        for (int i=0; i<n_nodes + 1; i++) {
            for (int a=0; a < n_actions; a++) {
                for (int j=0; j<n_nodes+1; j++) {
                    mdp.setTransitionProbability(i, a, j, 0.0);
                }
                //mdp.setTransitionProbability(i, a, i, 0.0);
            }
        }
        // no reward in the first state
        {
            Distribution* reward_density = 
                new SingularDistribution(0.0);
            for (int a=0; a<n_actions; a++) {
                mdp.setRewardDistribution(0,
                                          a,
                                          reward_density);
                }
        }

        for (int i=0; i<n_nodes; i++) {
            int n_edges = nodes[i]->outs.size();
            if (verbose >= 90) {
                printf ("Node %d has %d outgoing edges\n", i, n_edges);
            }
            for (int j=0; j<n_edges; j++) {
                Edge* edge = nodes[i]->outs[j];
                Distribution* reward_density = 
                    new SingularDistribution(edge->r);
                
                densities.push_back(reward_density);
                mdp.setTransitionProbability(i,
                                             edge->a,
                                             edge->dst->index,
                                             edge->p);
                if (verbose >= 90) {
                    printf ("%d: - a=%d - r=%f - p=%f ->%d\n", i, edge->a, edge->r, edge->p, edge->dst->index);
                }
                for (int a=0; a<n_actions; a++) {
                    mdp.setRewardDistribution(edge->dst->index,
                                              a,
                                              reward_density);
                }
            }
            
            if (!n_edges) {
                real mean_reward = nodes[i]->belief.getGreedyReturn(nodes[i]->state, gamma);
                for (int a=0; a<n_actions; a++) {
                    mdp.setTransitionProbability(i, a, terminal, 1.0);

                    Distribution* reward_density = 
                        new SingularDistribution(mean_reward);
                    densities.push_back(reward_density);
                    mdp.setRewardDistribution(i, a, reward_density);
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



/// A toy UCT stopping problem


//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int MakeDecision(int n_states, int n_actions, SimpleBelief prior, int state, real gamma, int n_iter, int verbose, int value_iterations, FILE* fout = NULL);

int main (int argc, char** argv)
{
    real alpha = 1.0;
    real beta = 1.0;
    int n_states = 2;
    int n_actions = 2;

    real gamma = 0.99;
    if (argc > 1) {
        gamma = atof(argv[1]);
    }

    int n_iter = 2;
    if (argc > 2) {
        n_iter = atoi(argv[2]);
    }

    int verbose = 0;
    if (argc > 3) {
        verbose = atoi(argv[3]);
    }

    int n_experiments = 1000;
    if (argc > 4) {
        n_experiments = atoi(argv[4]);
    }

    int horizon = 2.0/(1.0 - gamma);
    if (argc > 5) {
        horizon = atoi(argv[5]);
    }

    int value_iterations = n_iter + 1;
    if (argc > 6) {
        value_iterations = atoi(argv[6]);
    }

    

    std::vector<real> rewards(horizon);
    for (int t=0; t<horizon; t++) {
        rewards[t] = 0.0;
    }


    
    real average_oracle_return = 0.0;
    for (int experiment=0; experiment<n_experiments; experiment++) {
        real actual_probability = true_random(false); //urandom(0, 1);
        SimpleBelief belief(alpha, beta, -1.0, 1.0);
        int state = 0;
        for (int t=0; t<horizon; t++) {
            FILE* fout = NULL;
            if (n_experiments == 1) {
                char buffer[1024];
                sprintf(buffer, "test%d.dot", t);
                fout = fopen (buffer, "w");
                fprintf (fout, "digraph Lookahead {\n");
                fprintf (fout, "ranksep=2; rankdir=LR; \n");
            }

            int action = MakeDecision(n_states, n_actions, belief, state, gamma, n_iter, verbose, value_iterations, fout);

            if (fout) {
                fprintf (fout, "}\n");
                fclose (fout);
            }
            real reward = 0.0;
            
            if (state == 0 && action == 0) {
                if (urandom() < actual_probability) {
                    reward = 1.0;
                } else {
                reward = -1.0;
                }
            }
            
            rewards[t] += reward;
            int next_state = 0;
            if (action == 1) {
                next_state = 1;
            }
            
            belief.update(state, action, reward, next_state);
            state = next_state;
        }
        
        real oracle_return = 0.0;
        if (actual_probability > 0.5) {
            oracle_return = (2.0 * actual_probability - 1.0)/(1.0 - gamma);
        }
        average_oracle_return += oracle_return;
    }

    real inv_exp = 1.0 / (real) n_experiments;
    real total_reward = 0;
    real discount = 1.0;
    for (int t=0; t<horizon; t++) {
        rewards[t] *= inv_exp;
        total_reward += discount * rewards[t];
        discount *= gamma;
    }
    average_oracle_return *= inv_exp;
    

    
    printf("%f %d %f %f\n",
           gamma,
           n_iter,
           total_reward,
           average_oracle_return);
    return 0;
}

int MakeDecision(int n_states, int n_actions, SimpleBelief prior, int state, real gamma, int n_iter, int verbose, int value_iterations, FILE* fout)
{
    BeliefTree<SimpleBelief> tree(prior, state, n_states, n_actions);
    
    for (int iter=0; iter<n_iter; iter++) {
        std::vector<BeliefTree<SimpleBelief>::Node*> node_set = tree.getNodes();
        int edgeless_nodes = 0;
        int edge_nodes = 0;
        int node_index = -1;
        for (uint i=0; i<node_set.size(); ++i) {
            if (node_set[i]->outs.size()==0) {
                edgeless_nodes++;
                if (node_index == -1) {
                    node_index = i;
                }
            } else {
                edge_nodes++;
            }
        }

        if (verbose >= 100) {
            std::cout << edgeless_nodes << " leaf nodes, "
                      << edge_nodes << " non-leaf nodes, "
                      << "expanding node " << node_index
                      << std::endl;
        }

        if (node_index>=0) {
            tree.Expand(node_index, verbose); // VALGRIND
        } else {
            std::cout << "Warning: no nodes could be expanded\n";
        }
    }

    DiscreteMDP mdp = tree.CreateMDP(gamma, verbose); // VALGRIND
    mdp.Check();
    ValueIteration value_iteration(&mdp, gamma);
    if (verbose >= 75) {
        mdp.ShowModel();
    }

    
    double start_time = 0.0;
    double end_time = 0.0;

#if 1
    start_time = GetCPU();
    value_iteration.ComputeStateValues(0.00001, value_iterations);
    end_time = GetCPU();
    if (verbose) {
        std::cout << "# CPU " << end_time - start_time << std::endl;
    }

    if (verbose >= 60) {
        for (int s=0; s<mdp.GetNStates(); s++) {
            std::cout << "V[" << s << "]"
                      << " = " << value_iteration.getValue(s)
                      << std::endl;
        }
    }

    if (fout) {
        for (int s=0; s<mdp.GetNStates(); s++) {
            fprintf (fout,
                     "s%d [label = %2.f];\n",
                     s, value_iteration.getValue(s));
        }
    }
#endif

    start_time = GetCPU();
    value_iteration.ComputeStateActionValues(0.00001, value_iterations);
    end_time = GetCPU();
    if (verbose) {
        std::cout << "# CPU " << end_time - start_time << std::endl;
    }

    if (verbose >= 10) {
        for (int s=0; s<mdp.GetNStates(); s++) {
            for (int a=0; a<mdp.GetNActions(); a++) {
                std::cout << "Q[" << s << ", " << a << "]"
                          << " = " << value_iteration.getValue(s,a)
                      << std::endl;
            }
        }
    }

    if (fout) {
        mdp.dotModel(fout);
    }
    
    int a_max = 0;
    real Q_a_max = value_iteration.getValue(0, a_max);
    for (int a=1; a<n_actions; a++) {
        real Q_a = value_iteration.getValue(0, a);
        if (Q_a > Q_a_max) {
            a_max = a;
            Q_a_max = Q_a;
        }
    }

    return a_max;
}


#endif
