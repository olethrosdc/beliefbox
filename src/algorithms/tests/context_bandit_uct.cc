// -*- Mode: C++; -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/************************************************************************
 *                                                                      *
 * This program is free software; you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation; either version 2 of the License, or    *
 * (at your option) any later version.                                  *
 *                                                                      *
 ************************************************************************/

#ifdef MAKE_MAIN

#include "context_bandit_uct.h"

#include "UCB.h"
#include "Random.h"
#include "RandomNumberFile.h"
#include "MersenneTwister.h"
#include "ContextBandit.h"

/// A toy UCT stopping problem


//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int main (int argc, char** argv)
{
    //real alpha = 1.0;
    //real beta = 1.0;

    real mu_0 = 1.0;
    real tau_0 = 1.0;
    real tau = 0.5;
    //int n_states = 1;

    
    if (argc <=1) {
        fprintf (stderr, "Usage: bandit_uct method gamma n_iter n_actions verbose n_experiments horizon\n");
        exit(-1);
    }

    ExpansionMethod method = SerialExpansion;

    if (argc > 1) {
        method = (ExpansionMethod) atoi(argv[1]);
    }

    real gamma = 0.99;
    if (argc > 2) {
        gamma = atof(argv[2]);
    }

    int n_iter = 2;
    if (argc > 3) {
        n_iter = atoi(argv[3]);
    }

    int n_actions = 2;
    if (argc > 4) {
        n_actions = atoi(argv[4]);
        assert(n_actions > 1);
    }

    int verbose = 0;
    if (argc > 5) {
        verbose = atoi(argv[5]);
    }

    int n_experiments = 1000;
    if (argc > 6) {
        n_experiments = atoi(argv[6]);
    }

    int horizon = (int) floor(2.0/(1.0 - gamma));
    if (argc > 7) {
        horizon = atoi(argv[7]);
    }

    int max_value_iterations = n_iter + 1;

    std::vector<real> rewards(horizon);
    for (int t=0; t<horizon; t++) {
        rewards[t] = 0.0;
    }
    
    real average_oracle_return = 0.0;
    real average_oracle_reward = 0.0;

    //RandomNumberFile rng("./dat/r1e7.bin");
    int n_states  = 1;

    srand(1228517343);
    srand48(1228517343);
    
    MersenneTwisterRNG rng;
    rng.manualSeed(18517339);
    // perform experiments
    real start_time = GetCPU();

        //Create an empty tree
    ContextBanditBelief prior_belief(n_states, n_actions, tau, mu_0, tau_0);
    

    for (int experiment=0; experiment<n_experiments; experiment++) {
        BeliefTree<ContextBanditBelief>* new_tree= new BeliefTree<ContextBanditBelief> (prior_belief, 0, n_states, n_actions, gamma);

        ContextBandit context_bandit(n_actions, 3, 8, &rng);

        // initial state and belief
        int state = 0;
        ContextBanditBelief belief(n_states, n_actions, tau, mu_0, tau_0);
        std::vector<real> Er(2);

        real oracle_return = 0.0;

        // looop over time
        for (int t=0; t<horizon; t++) {
            // write graph to a file if we are only doing one experiment
            FILE* fout = NULL;
            if (n_experiments == 1 && verbose > 50) {
                char buffer[1024];
                sprintf(buffer, "test%d.dot", t);
                fout = fopen (buffer, "w");
                fprintf (fout, "digraph Lookahead {\n");
                fprintf (fout, "ranksep=2; rankdir=LR; \n");
            }
            
            state = context_bandit.getState();

            int action = MakeDecision(*new_tree,
                                      method,
                                      n_states,
                                      n_actions,
                                      belief,
                                      state,
                                      gamma,
                                      n_iter,
                                      verbose,
                                      max_value_iterations,
                                      fout);

            // close file
            if (fout) {
                fprintf (fout, "}\n");
                fclose (fout);
            }
            
            for (int i=0; i<n_actions; i++) {
                Er[i] = context_bandit.getMean(i);
            }
            oracle_return += Max(Er);
            
            // act
            context_bandit.Act(action);

            // calculate reward
            real reward = context_bandit.getReward();
            rewards[t] += reward;

            // calculate next state
            int next_state = context_bandit.getState(); // no change in state

            
            if ((n_experiments == 1) && (verbose >= 90)) {
                printf("#SARS: %d %d -> %f %d\n",
                       state, action, reward, next_state);
            }

            // update belief
            belief.update(state, action, reward, next_state);

            // update future belief tree
            BeliefTreeNode* new_root = new_tree->FindObservation(new_tree->root,
                                                                 action,
                                                                 reward,
                                                                 next_state);
            
            if (new_root) {
                new_tree->MakeRoot(new_root, verbose);
            } else {
                delete new_tree;
                new_tree = new BeliefTree<ContextBanditBelief> (belief, next_state, n_states, n_actions, gamma);
            }

            // update state
            state = next_state;
        } // for t
        
        delete new_tree;

        // collect statistics

        average_oracle_return += oracle_return;
        average_oracle_reward += ((real) horizon);

    } // for experiment

    // average statistics
    real inv_exp = 1.0 / (real) n_experiments;
    real total_reward = 0.0;
    real discounted_reward = 0.0;
    real discount = 1.0;
    for (int t=0; t<horizon; t++) {
        //rewards[t] *= inv_exp;
        discounted_reward += discount * rewards[t];
        total_reward += rewards[t];
        discount *= gamma;
    }
    average_oracle_return *= inv_exp;
    average_oracle_reward *= inv_exp;
    discounted_reward *= inv_exp;
    total_reward *= inv_exp;

    real end_time = GetCPU();
    // show statistics
    printf("%d %d %f %d %f %f %f %f %f %d %d\n",
           method,
           n_iter,
           gamma,
           n_actions,
           total_reward,
           discounted_reward,
           average_oracle_reward,
           average_oracle_return,
           end_time - start_time,
           n_experiments,
           horizon);

    FILE* fstats = fopen("reward.out", "w");
    for (int t=0; t<horizon; t++) {
        fprintf(fstats, "%f\n",rewards[t]*inv_exp);
    }
    fclose(fstats);
    return 0;
}



                    
                    
int ChooseGreedyBranch(BeliefTree<ContextBanditBelief>& tree, real gamma)
{
    DiscreteMDP mean_mdp = tree.CreateMeanMDP(gamma, 0);
    mean_mdp.Check();

    DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, 0);
    upper_mdp.Check();

    ValueIteration value_iteration_lower(&mean_mdp, gamma);
    ValueIteration value_iteration_upper(&upper_mdp, gamma);
    return 0;
    
}

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
                 FILE* fout)
{

    double start_time = GetCPU();

    double vi_threshold = 0.0001;
    //BeliefTree<ContextBanditBelief> tree(prior, state, n_states, n_actions, gamma);
    BeliefTree<ContextBanditBelief> &tree = new_tree;

    int n_rewards = 2;
    int branching_factor = n_actions * n_states * n_rewards;
    int starting_iter = (tree.GetNumberOfLeafNodes() - 1)/(branching_factor - 1);

    for (int iter = starting_iter; iter<n_iter; iter++) {
        BTNodeSet node_set = tree.getNodes();
        int n_edge_nodes = 0;
        BeliefTreeNode* selected_node = NULL;
        std::vector<BeliefTreeNode*> leaf_nodes;
        
        // leaf ndoes are state nodes
        for (BTNodeSet::iterator i=node_set.begin(); i!=node_set.end(); ++i) {
            BeliefTreeNode* node = *i;
            if (node->outs.size() == 0) {
                leaf_nodes.push_back(node);
            } else {
                n_edge_nodes++;
            }
        }
        
        int n_leaf_nodes = (int) leaf_nodes.size();
        
        // Find the root actions
        std::vector<int> root_action(n_leaf_nodes);
        for (int i=0; i<n_leaf_nodes; i++) {
            root_action[i] = tree.FindRootAction(leaf_nodes[i]);
        }
        
        if (expansion_method == SerialExpansion) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                std::vector<real> &Ub = node->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node->L = Li;
                if (Li > Ui) {
                    Ui = Li;
                }
                U[i] = Ui;
            }
            selected_node = leaf_nodes[0];
        } else if (expansion_method == RandomExpansion) {
            // randomly expand the tree
            int X =  (int) floor(urandom()*((real) leaf_nodes.size()));
            selected_node = leaf_nodes[X];
        } else if (expansion_method == HighestMeanValue) {
            // return the node with the highest lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node->L = Li;
                U[i] = Li;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == HighestDiscountedMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                U[i] = node->R + pow(gamma, (real) node->depth)*Li;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                U[i] = node->belief.sampleReturn(node->state, gamma);
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonBound) {
            std::vector<real> U(leaf_nodes.size());
            //std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] = Ui;
                //L[i] = Li;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        }  else if (expansion_method == DiscountedThompsonBound) {
            std::vector<real> U(leaf_nodes.size());
            //std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  pow(gamma, (real) node->depth) * Ui;
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
                //L[i] = Li;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == HighProbabilityBoundOnly) {
            // returns the maximum of the current upper bounds
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                node->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node->U);
                U[i] = Ui;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedHighProbabilityBoundOnly) {
            // returns the maximum of the current upper bounds
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                node->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node->U);
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
            }
            selected_node = leaf_nodes[ArgMax(U)]; 
        } else if (expansion_method == HighProbabilityBound) {
            // returns the maximum of the current upper bounds or the lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                node->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node->U);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  Ui;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedHighProbabilityBound) {
            // returns the maximum of the current upper bounds or the lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                node->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node->U);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
                //U[i] = pow(gamma, (real) node->depth) * Ui;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == MeanHighProbabilityBound) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                std::vector<real> &Ub = node->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node->L = Li;
                U[i] = Ui;
            }
            selected_node = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedMeanHighProbabilityBound) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                std::vector<real> &Ub = node->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Li > Ui) {
                    Ui = Li;
                }
                //real p = node->GetPathProbability();
                U[i] = (node->R + pow(gamma, (real) node->depth * Ui));
            }
            selected_node = leaf_nodes[ArgMax(U)]; 
        } else if (expansion_method == GreedyBoundReduction) {
            // Find the node with the highest bound
            std::vector<real> U(leaf_nodes.size());
            std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                std::vector<real> &Ub = node->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real discount = node->GetPathProbability() * pow(gamma, (real) node->depth);
                real Ui = discount*Mean(Ub);
                real Li = discount*node->belief.getGreedyReturn(node->state, gamma);

                //
                //if (Li > Ui) {
                //Ui = Li;
                //}
                
                U[i] = Ui;
                L[i] = Li;
            }

            selected_node = leaf_nodes[ArgMax(U)];


            int argmax_L = ArgMax(L);
            real max_L = L[argmax_L];
            selected_node = leaf_nodes[argmax_L];
            int best_action = root_action[argmax_L];


            // Find the highest upper bound
            int argmax_U = -1;
            real max_U = -INF;
            for (uint i=0; i<U.size(); i++) {
                if (root_action[i] != best_action) {
                    if (L[i] < max_L) {
                        if ((argmax_U == -1) || (max_U < U[i])) {
                            max_U = U[i];
                            argmax_U = i;
                        }
                    }
                }
            }

            // if highest upper bound beats highest lower bound
            // then try to expand the highest upper bound
            // otherwise try to expand the lower bound..
            if (max_U > max_L) {
                selected_node = leaf_nodes[argmax_U];
            } else {
                selected_node = leaf_nodes[argmax_L];
            }

        } else if (expansion_method == BAST) {
                // sample all leaf nodes
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                //assert(node->index==n);
                std::vector<real> &Ub = node->U;
                Ub.push_back(node->belief.sampleReturn(node->state, gamma));
            }
                // propagate upper bounds to the root
            DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, verbose);
            upper_mdp.Check(); 
            ValueIteration value_iteration_upper(&upper_mdp, gamma);
            value_iteration_upper.ComputeStateActionValues(vi_threshold, max_value_iterations);

            
            int s = 0; // the MDP state
            while (1) {
                std::vector<real> Ua(n_actions);
                for (int a=0; a<n_actions; a++) {
                    Ua[a] = value_iteration_upper.getValue(s,a);
                }
                int a_max = ArgMax(Ua);
                // TODO Choose only ONE node via sampling (or multiple nodes?)
                // Warning: uses STATE value as ID
                BeliefTreeNode* node = NULL;
                for (BTNodeSet::iterator i=node_set.begin(); i!=node_set.end(); ++i) {
                    BeliefTreeNode* inode = *i;
                    if (inode->index == s) {
                        node = inode;
                    }
                }
                if (!node) {
                    fprintf(stderr, "Error: Could not find node!\n");
                    exit(-1);
                }

                //printf ("state: %d selected_node: %d with %d edges\n", s, node->index, node->outs.size());
                // if node is a leaf node...
                if (node->outs.size() == 0) {
                    selected_node = node;
                    break;
                }

                // if there are edges, randomly select one of the edges
                // to traverse
                //std::vector<BeliefTreeEdge*> edges;
                real pr = 0.0;
                real X = urandom();
                for (BTEdgeSet::iterator j = node->outs.begin();
                     j != node->outs.end(); ++j) {
                    BeliefTreeEdge* edge = *j;
                    if (edge->a == a_max) {
                        //edges.push_back(edge);
                        pr += edge->p;
                        if (pr >= X) {
                            s = edge->dst->index;
                            break;
                        }
                    }
                }
                //printf("Total probability :%f\n", pr);
            }
            
            //printf("exiting\n");
       } else if (expansion_method == BAST_FULL) {
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = leaf_nodes[i];
                //assert(node->index==n);
                std::vector<real> &Ub = node->U;
                Ub.push_back(node->belief.sampleReturn(node->state, gamma));
            }
                // propagate upper bounds to the root
            DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, verbose);
            upper_mdp.Check(); 
            ValueIteration value_iteration_upper(&upper_mdp, gamma);
            value_iteration_upper.ComputeStateActionValues(vi_threshold, max_value_iterations);

            
            int s = 0; // the MDP state
            while (1) {
                std::vector<real> Ua(n_actions);
                for (int a=0; a<n_actions; a++) {
                    Ua[a] = value_iteration_upper.getValue(s,a);
                }
                int a_max = ArgMax(Ua);
                // TODO Choose only ONE node via sampling
                // (or multiple nodes?)
                // Warning: uses STATE value as ID
                BeliefTreeNode* node = NULL;
                for (BTNodeSet::iterator i=node_set.begin(); i!=node_set.end(); ++i) {
                    BeliefTreeNode* inode = *i;
                    if (inode->index == s) {
                        node = inode;
                    }
                }
                if (!node) {
                    fprintf(stderr, "Error: Could not find node!\n");
                    exit(-1);
                }

                //printf ("state: %d selected_node: %d with %d edges\n", s, node->index, node->outs.size());
                // if node is a leaf node...
                if (node->outs.size() == 0) {
                    selected_node = node;
                    break;
                }

                // if there are edges, randomly select one of the edges
                // to traverse
                //std::vector<BeliefTreeEdge*> edges;
                real pr = 0.0;
                real X = urandom();
                for (BTEdgeSet::iterator j = node->outs.begin();
                     j != node->outs.end(); ++j) {
                    BeliefTreeEdge* edge = *j;
                    if (edge->a == a_max) {
                        //edges.push_back(edge);
                        pr += edge->p;
                        if (pr >= X) {
                            s = edge->dst->index;
                            break;
                        }
                    }
                }
                //printf("Total probability :%f\n", pr);
            }
            
            //printf("exiting\n");
        } else {
            std::cerr << "Unknown method " << expansion_method << std::endl;
            exit(-1);
        }

        if (verbose >= 100) {
            std::cout << n_leaf_nodes << " leaf nodes, "
                      << n_edge_nodes << " edge nodes, "
                      << node_set.size() << " total nodes, "
                      << "expanding node " << selected_node->index
                      << std::endl;
        }

        if (selected_node) {
            tree.Expand(selected_node, verbose);
        } else {
            std::cout << "Warning: no nodes could be expanded\n";
        }

        //std::cout << node_set.size() << " total nodes"
        //              << std::endl;
    } // for(iter)

    if (verbose >= 90) {
        std::cout << "Final expansion : "
                  << tree.GetNumberOfLeafNodes() << " leaf nodes, " 
                  << tree.getNodes().size() << " nodes"
                  << std::endl;
    }
    BTNodeSet& node_set = tree.getNodes();
    for (BTNodeSet::iterator i=node_set.begin(); i!=node_set.end(); ++i) {
        BeliefTreeNode* node = *i;
        node->L = node->belief.getGreedyReturn(node->state, gamma);
    }

    DiscreteMDP mean_mdp = tree.CreateMeanMDP(gamma, verbose);
    mean_mdp.Check();

    DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, verbose);
    upper_mdp.Check();

    ValueIteration value_iteration_lower(&mean_mdp, gamma);
    ValueIteration value_iteration_upper(&upper_mdp, gamma);
    if (verbose >= 75) {
        mean_mdp.ShowModel();
        upper_mdp.ShowModel();
    }


#if 1
    if (fout) {
        value_iteration_lower.ComputeStateValues(vi_threshold, max_value_iterations);
        value_iteration_upper.ComputeStateValues(vi_threshold, max_value_iterations);
        
        if (verbose >= 60) {
            std::cout << "mdp states " << mean_mdp.getNStates()
                      << " nodes " << node_set.size()
                      << std::endl;
            for (int s=0; s<mean_mdp.getNStates() - 1; s++) {
                if (s >= (int) node_set.size()) {
                    std::cerr << "ERROR in mdp size!"
                              << s << " >= " << node_set.size()
                              << "mdp size = " << mean_mdp.getNStates()
                              << std::endl;
                    exit(-1);
                }
                real L_bound = -INF;
                real U_bound = INF;
                BeliefTreeNode* node = NULL;
                for (BTNodeSet::iterator i=node_set.begin(); i!=node_set.end(); ++i) {
                    BeliefTreeNode* inode = *i;
                    if (inode->index == s) {
                        node = inode;
                    }
                }  

                L_bound = node->belief.getGreedyReturn(node->state, gamma);
                U_bound = Mean(node->U);
                real d = (real) (1+node->depth);
                real gd = pow(gamma, d);
                if (isnan(U_bound)) {
                    U_bound = INF;
                }
                L_bound *= gd;
                U_bound *= gd;
                if (L_bound > value_iteration_lower.getValue(s)) {
                    std::cout << "+";
                }
                if (U_bound < value_iteration_upper.getValue(s)) {
                    std::cout << "-";
                }
                
                std::cout << "V(" << s << ")"
                          << " in [" << value_iteration_lower.getValue(s)
                          << ", " << value_iteration_upper.getValue(s)
                          << "] C [ " << L_bound
                          << ", " << U_bound
                          << "] " << std::endl;
            }
        }


        for (int s=0; s<mean_mdp.getNStates(); s++) {
            fprintf (fout,
                     "s%d [label = \"%.2f - %.2f\"];\n",
                     s,
                     value_iteration_lower.getValue(s),
                     value_iteration_upper.getValue(s));
        }
    }
#endif


    value_iteration_upper.ComputeStateActionValues(vi_threshold, max_value_iterations);
    value_iteration_lower.ComputeStateActionValues(vi_threshold, max_value_iterations);
    double end_time = GetCPU();
    if (verbose) {
        std::cout << "# CPU " << end_time - start_time << std::endl;
    }

    if (verbose >= 20) {
        for (int s=0; s<mean_mdp.getNStates(); s++) {
            for (int a=0; a<mean_mdp.getNActions(); a++) {
                std::cout << "Q[" << s << ", " << a << "]"
                          << " = " << value_iteration_lower.getValue(s,a)
                          << " - " << value_iteration_upper.getValue(s,a)
                          << std::endl;
            }
        }
    }





    if (fout) {
        mean_mdp.dotModel(fout);
    }
    
    int a_maxL = 0;
    real Q_a_maxL = value_iteration_lower.getValue(0, a_maxL);
    for (int a=1; a<n_actions; a++) {
        real Q_a = value_iteration_lower.getValue(0, a);
        if (Q_a > Q_a_maxL) {
            a_maxL = a;
            Q_a_maxL = Q_a;
        }
    }

#if 1
    int a_maxU = -1;
    real Q_a_maxU = -INF;
    std::vector<real> VU(n_actions);
    for (int a=0; a<n_actions; a++) {
        real Q_a = value_iteration_upper.getValue(0, a);
        VU[a] = Q_a;
        if ((a != a_maxL)
            &&
            ((Q_a > Q_a_maxU) || (a==0))
            ) {
            a_maxU = a;
            Q_a_maxU = Q_a;
        }
    }

    if (verbose >= 5) {
        std::cout << Q_a_maxL - Q_a_maxU << " "
                  << Q_a_maxU - value_iteration_lower.getValue(a_maxU) << " "
                  << value_iteration_upper.getValue(a_maxL) - Q_a_maxL << " "
                  << " # DQ" << std::endl;
    }
#endif    
    return a_maxL;
    //return ArgMax(VU);
}



#endif
