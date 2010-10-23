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

#include "bandit_uct.h"

#include "UCB.h"
#include "Random.h"
#include "RandomNumberFile.h"

/// A toy UCT stopping problem

/* void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);
   TODO: implement algorithm
*/

int main (int argc, char** argv)
{
    real alpha = 1.0;
    real beta = 1.0;
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
    int n_states  = 0;

    srand48(1228517343);
    // perform experiments
    real start_time = GetCPU();
    for (int experiment=0; experiment<n_experiments; experiment++) {
        
        std::vector<real> Er(n_actions);

        for (int i=0; i<n_actions; i++) {
            Er[i] = drand48();//rng.uniform();
            //printf ("%f ", Er[i]);
        }
        //printf("#ER[] \n");

        // initial state and belief
        int state = 0;
        BanditBelief belief(n_actions, alpha, beta);
                                                
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

            if (verbose >= 10) {
                printf(" Er:");
                for (int i=0; i<n_actions; i++) {
                    printf(" %f", Er[i]);
                }
                printf("\nhEr:");
                for (int i=0; i<n_actions; i++) {
                    printf(" %f", belief.getMeanReward(0, i));
                }
                printf("\n");
            }
            int action = MakeDecision(method,
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

            // calculate reward
            real reward = 0.0;
            if (urandom() < Er[action]) {
                reward = 1.0;
            } else {
                reward = 0.0;

            }
            rewards[t] += reward;

            // calculate next state
            int next_state = 0; // no change in state

            
            if ((n_experiments == 1) && (verbose >= 90)) {
                printf("#SARS: %d %d -> %f %d\n",
                       state, action, reward, next_state);
            }
            // update belief
            belief.update(state, action, reward, next_state);

            // update state
            state = next_state;
        }
        
        // collect statistics
        real oracle_return = 0.0;
        real max_reward = Max(Er);
        oracle_return = max_reward*(1.0 - pow(gamma, (real) (horizon+1)))/(1.0 - gamma);
        average_oracle_return += oracle_return;
        average_oracle_reward += max_reward * ((real) horizon);

    }

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
    printf("%d %d %f %d %f %f %f %f %f %d\n",
           method,
           n_iter,
           gamma,
           n_actions,
           total_reward,
           discounted_reward,
           average_oracle_reward,
           average_oracle_return,
           end_time - start_time,
           n_experiments);
    return 0;
}

#if 0
static bool sort_using_greater_than(real u, real v)
{
    return u > v;
}
#endif


                    
                    
int ChooseGreedyBranch(BeliefTree<BanditBelief>& tree, real gamma)
{
    DiscreteMDP mean_mdp = tree.CreateMeanMDP(gamma, 0);
    mean_mdp.Check();

    DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, 0);
    upper_mdp.Check();

    ValueIteration value_iteration_lower(&mean_mdp, gamma);
    ValueIteration value_iteration_upper(&upper_mdp, gamma);
    return 0;
    
}

int MakeDecision(ExpansionMethod expansion_method,
                 int n_states,
                 int n_actions,
                 BanditBelief prior,
                 int state,
                 real gamma,
                 int n_iter,
                 int verbose,
                 int max_value_iterations,
                 FILE* fout)
{

    double start_time = GetCPU();

    BeliefTree<BanditBelief> tree(prior, state, n_states, n_actions, gamma);


    for (int iter=0; iter<n_iter; iter++) {
        std::vector<BeliefTreeNode*> node_set = tree.getNodes();
        int n_edge_nodes = 0;
        int node_index = -1;
        int action_node_index = -1;
        std::vector<int> leaf_nodes;
        std::set<int> action_leaf_nodes;
        
        // leaf ndoes are state nodes
        for (uint i=0; i<node_set.size(); ++i) {
            if (node_set[i]->outs.size()==0) {
                leaf_nodes.push_back(i);
            } else {
                n_edge_nodes++;
            }
        }
        
        // do this when leaf nodes are action nodes
        for (uint i=0; i<leaf_nodes.size(); ++i) {
            BeliefTreeNode* parent = node_set[leaf_nodes[i]]->GetParent();
            if (parent) {
                for (uint j=0; j<node_set.size(); ++j) {
                    if (parent == node_set[j]) {
                        action_leaf_nodes.insert(j);
                        break;
                    } 
                }
            } else {
                action_leaf_nodes.insert(leaf_nodes[i]);
            }
        }

        int n_leaf_nodes = (int) leaf_nodes.size();
        
        // Find the root actions
        std::vector<int> root_action(n_leaf_nodes);
        for (int i=0; i<n_leaf_nodes; i++) {
            root_action[i] = tree.FindRootAction(node_set[leaf_nodes[i]]);
        }

        
        if (expansion_method == SerialExpansion) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                std::vector<real> &Ub = node_set[n]->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node_set[n]->L = Li;
                if (Li > Ui) {
                    Ui = Li;
                }
                U[i] = Ui;
            }
            node_index = leaf_nodes[0];
            action_node_index = *action_leaf_nodes.begin();
        } else if (expansion_method == RandomExpansion) {
            int X =  (int) floor(urandom()*((real) leaf_nodes.size()));
            node_index = leaf_nodes[X];
#if 0

            
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 ++i, --X) {
                if (X==0) {
                    action_node_index = *i;
                    break;
                }
            }
#endif
            //node_index = leaf_nodes[X];
        } else if (expansion_method == HighestMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                std::vector<real> &Ub = node_set[n]->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node_set[n]->L = Li;
                if (Li > Ui) {
                    Ui = Li;
                }
                U[i] = Li;
            }

            node_index = leaf_nodes[ArgMax(U)];
#if 0
            DiscreteMDP mean_mdp = tree.CreateMeanMDP(gamma, verbose);
            mean_mdp.Check();
            ValueIteration value_iteration_lower(&mean_mdp, gamma);
            value_iteration_lower.ComputeStateValues(0.00001, max_value_iterations);
            std::vector<real> Ua(action_leaf_nodes.size());
            int j=0;
            real Ua_max = 0.0;
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 ++i, ++j) {
                int index = *i;
                BeliefTreeNode* node = node_set[index];
                Ua[j] = pow(gamma, (real) node->depth) *value_iteration_lower.getValue(index);
                if ((action_node_index == -1) || (Ua_max < Ua[j])) {
                    Ua_max = Ua[j];
                    action_node_index = index;
                }
            }
            //action_node_index = *action_leaf_nodes.begin();
            std::cout << "action indices: ";
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 i++) {
                std::cout << *i << " ";
            }
            std::cout << "\nselected index: " << action_node_index
                      << std::endl;
#endif
        } else if (expansion_method == HighestDiscountedMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                std::vector<real> &Ub = node_set[n]->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node_set[n]->L = Li;
                if (Li > Ui) {
                    Ui = Li;
                }
                U[i] = node->R + pow(gamma, (real) node->depth)*Li;
            }
            node_index = leaf_nodes[ArgMax(U)];
#if 0
            DiscreteMDP mean_mdp = tree.CreateMeanMDP(gamma, verbose);
            mean_mdp.Check();
            ValueIteration value_iteration_lower(&mean_mdp, gamma);
            value_iteration_lower.ComputeStateValues(0.00001, max_value_iterations);
            std::vector<real> Ua(action_leaf_nodes.size());
            int j=0;
            real Ua_max = 0.0;
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 ++i, ++j) {
                int index = *i;
                BeliefTreeNode* node = node_set[index];
                Ua[j] = pow(gamma, (real) node->depth) *value_iteration_lower.getValue(index);
                if ((action_node_index == -1) || (Ua_max < Ua[j])) {
                    Ua_max = Ua[j];
                    action_node_index = index;
                }
            }

            //action_node_index = *action_leaf_nodes.begin();
            std::cout << "action indices: ";
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 i++) {
                std::cout << *i << " ";
            }
            std::cout << "\nselected index: " << action_node_index
                      << std::endl;
#endif
        } else if (expansion_method == ThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = node_set[leaf_nodes[i]];
                U[i] = node->belief.sampleReturn(node->state, gamma);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = node_set[leaf_nodes[i]];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonBound) {
            std::vector<real> U(leaf_nodes.size());
            //std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = node_set[leaf_nodes[i]];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] = Ui;
                //L[i] = Li;
            }
            node_index = leaf_nodes[ArgMax(U)];
        }  else if (expansion_method == DiscountedThompsonBound) {
            std::vector<real> U(leaf_nodes.size());
            //std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTreeNode* node = node_set[leaf_nodes[i]];
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
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == HighProbabilityBoundOnly) {
            // returns the maximum of the current upper bounds
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                U[i] = Ui;
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedHighProbabilityBoundOnly) {
            // returns the maximum of the current upper bounds
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
            }
            node_index = leaf_nodes[ArgMax(U)]; 
        } else if (expansion_method == HighProbabilityBound) {
            // returns the maximum of the current upper bounds or the lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  Ui;
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedHighProbabilityBound) {
            // returns the maximum of the current upper bounds or the lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                real p = node->GetPathProbability();
                U[i] = p*(node->R + pow(gamma, (real) node->depth) * Ui);
                //U[i] = pow(gamma, (real) node->depth) * Ui;
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == MeanHighProbabilityBound) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                std::vector<real> &Ub = node->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                node_set[n]->L = Li;
                // actually makes no difference.
                //if (Li > Ui) {
                //Ui = Li;
                //}
                U[i] = Ui;
            }
            node_index = leaf_nodes[ArgMax(U)];
#if 0
            DiscreteMDP mean_mdp = tree.CreateMeanMDP(gamma, verbose);
            mean_mdp.Check();
            ValueIteration value_iteration_lower(&mean_mdp, gamma);
            value_iteration_lower.ComputeStateValues(0.00001, max_value_iterations);

            DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, verbose);
            upper_mdp.Check();
            ValueIteration value_iteration_upper(&upper_mdp, gamma);
            value_iteration_upper.ComputeStateValues(0.00001, max_value_iterations);
            std::vector<real> Ua(action_leaf_nodes.size());
            int j=0;
            real Ua_max = 0.0;
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 ++i, ++j) {
                int index = *i;
                //BeliefTreeNode* node = node_set[index];
                Ua[j] = value_iteration_upper.getValue(index);
                if ((action_node_index == -1) || (Ua_max < Ua[j])) {
                    Ua_max = Ua[j];
                    action_node_index = index;
                }
            }
#endif
        } else if (expansion_method == DiscountedMeanHighProbabilityBound) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
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
            node_index = leaf_nodes[ArgMax(U)]; 
#if 0
            DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, verbose);
            upper_mdp.Check();
            ValueIteration value_iteration_upper(&upper_mdp, gamma);
            value_iteration_upper.ComputeStateValues(0.00001, max_value_iterations);
            std::vector<real> Ua(action_leaf_nodes.size());
            int j=0;
            real Ua_max = 0.0;
            for (std::set<int>::iterator i=action_leaf_nodes.begin();
                 i!=action_leaf_nodes.end();
                 ++i, ++j) {
                int index = *i;
                BeliefTreeNode* node = node_set[index];
                Ua[j] = node->R + pow(gamma, (real) node->depth) * value_iteration_upper.getValue(index);
                if ((action_node_index == -1) || (Ua_max < Ua[j])) {
                    Ua_max = Ua[j];
                    action_node_index = index;
                }
            }
#endif
        } else if (expansion_method == GreedyBoundReduction) {
            // Find the node with the highest bound
            std::vector<real> U(leaf_nodes.size());
            std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
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

            node_index = leaf_nodes[ArgMax(U)];


            int argmax_L = ArgMax(L);
            real max_L = L[argmax_L];
            node_index = leaf_nodes[argmax_L];
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
                node_index = leaf_nodes[argmax_U];
            } else {
                node_index = leaf_nodes[argmax_L];
            }
        } else if (expansion_method == BAST) {
                // sample all leaf nodes
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTreeNode* node = node_set[n];
                assert(node->index==n);
                std::vector<real> &Ub = node_set[n]->U;
                Ub.push_back(node->belief.sampleReturn(node->state, gamma));
            }
                // propagate upper bounds to the root
            DiscreteMDP upper_mdp = tree.CreateUpperBoundMDP(gamma, verbose);
            upper_mdp.Check(); 
            ValueIteration value_iteration_upper(&upper_mdp, gamma);
            value_iteration_upper.ComputeStateActionValues(0.00001, max_value_iterations);

            
            int s = 0;
            while (1) {
                std::vector<real> Ua(n_actions);
                for (int a=0; a<n_actions; a++) {
                    Ua[a] = value_iteration_upper.getValue(s,a);
                }
                int a_max = ArgMax(Ua);
                // TODO: Choose only ONE node via sampling (or multiple nodes?)
                BeliefTreeNode* node = node_set[s];
                    //printf ("state: %d node_index: %d with %d edges\n", s, node->index, node->outs.size());
                if (node->outs.size() == 0) {
                    node_index = node->index;
                    break;
                }

                    // if there are edges, randomly select one of the edges
                    // to traverse
                std::vector<BeliefTreeEdge*> edges;
                real pr = 0.0;
                real X = urandom();
                for (uint i=0; i<node->outs.size(); i++) {
                    BeliefTreeEdge* edge = node->outs[i];
                        //printf("node %d edge %d: %d %f\n", node->index, i, edge->a, edge->p);
                    if (edge->a == a_max) {
                        edges.push_back(edge);
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
                      << "expanding node " << node_index
                      << " & action node " << action_node_index
                      << std::endl;
        }

        if (node_index>=0) {
            tree.Expand(node_index, verbose);
        } else if (action_node_index>=0) {
            BeliefTreeNode* action_node = node_set[action_node_index];
            //std::cout << action_node->outs.size() << " descendants\n";
            if (action_node->outs.size() == 0) {
                tree.Expand(action_node_index, verbose);
            }
            for (uint e=0; e<action_node->outs.size(); ++e) {
                BeliefTreeNode* dst = action_node->outs[e]->dst;
                
                for (uint d=0; d<leaf_nodes.size(); ++d) {
                    if (node_set[leaf_nodes[d]] == dst) {
                        assert(dst->outs.size()==0);
                        if (verbose >= 100) {
                            std::cout << "-> descendant: " << leaf_nodes[d]
                                      << std::endl;
                        }
                        tree.Expand(leaf_nodes[d], verbose);
                    }
                }
            }
        } else {
            std::cout << "Warning: no nodes could be expanded\n";
        }

        //std::cout << node_set.size() << " total nodes"
        //              << std::endl;
    } // for(iter)


    std::vector<BeliefTreeNode*> node_set = tree.getNodes();
    for (uint s=0; s<node_set.size(); s++) {
        node_set[s]->L = node_set[s]->belief.getGreedyReturn(node_set[s]->state, gamma);
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
        value_iteration_lower.ComputeStateValues(0.00001, max_value_iterations);
        value_iteration_upper.ComputeStateValues(0.00001, max_value_iterations);
        
        if (verbose >= 60) {
            std::cout << "mdp states " << mean_mdp.GetNStates()
                      << " nodes " << node_set.size()
                      << std::endl;
            for (int s=0; s<mean_mdp.GetNStates(); s++) {
                real L_bound = -INF;
                real U_bound = INF;
                if (s < (int) node_set.size()) {
                    //L_bound = node_set[s]->L;
                    L_bound = node_set[s]->belief.getGreedyReturn(node_set[s]->state, gamma);
                    U_bound = Mean(node_set[s]->U);
                    real d = (real) (1+node_set[s]->depth);
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
                }
                std::cout << "V(" << s << ")"
                          << " in [" << value_iteration_lower.getValue(s)
                          << ", " << value_iteration_upper.getValue(s)
                          << "] C [ " << L_bound
                          << ", " << U_bound
                          << "] " << std::endl;
            }
        }

        for (int s=0; s<mean_mdp.GetNStates(); s++) {
            fprintf (fout,
                     "s%d [label = \"%.2f - %.2f\"];\n",
                     s,
                     value_iteration_lower.getValue(s),
                     value_iteration_upper.getValue(s));
        }
    }
#endif


    value_iteration_upper.ComputeStateActionValues(0.00001, max_value_iterations);
    value_iteration_lower.ComputeStateActionValues(0.00001, max_value_iterations);
    double end_time = GetCPU();
    if (verbose) {
        std::cout << "# CPU " << end_time - start_time << std::endl;
    }

    if (verbose >= 20) {
        for (int s=0; s<mean_mdp.GetNStates(); s++) {
            for (int a=0; a<mean_mdp.GetNActions(); a++) {
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

#if 0
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
        std::cout << "DQ = " << Q_a_maxL - Q_a_maxU << std::endl;
    }
#endif    
    return a_maxL;
    //return ArgMax(VU);
}



#endif
