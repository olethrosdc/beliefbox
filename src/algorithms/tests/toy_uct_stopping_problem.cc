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

#include "toy_uct_stopping_problem.h"
#include "RandomNumberFile.h"

/// A toy UCT stopping problem


//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int main (int argc, char** argv)
{
    real alpha = 1.0;
    real beta = 1.0;
    int n_states = 2;
    int n_actions = 2;


    enum ExpansionMethod expansion_method = SerialExpansion;
    if (argc > 1) {
        expansion_method = (enum ExpansionMethod) atoi(argv[1]);
    } else {
        printf ("Usage: toy_uct_stopping_problem method gamma n_iter verbose experiments horizon max_value_iterations\n");
    }
    
    real gamma = 0.99;
    if (argc > 2) {
        gamma = atof(argv[2]);
    }

    int n_iter = 2;
    if (argc > 3) {
        n_iter = atoi(argv[3]);
    }

    int verbose = 0;
    if (argc > 4) {
        verbose = atoi(argv[4]);
    }

    int n_experiments = 1000;
    if (argc > 5) {
        n_experiments = atoi(argv[5]);
    }

    int horizon = 2.0/(1.0 - gamma);
    if (argc > 6) {
        horizon = atoi(argv[6]);
    }

    int max_value_iterations = n_iter + 1;
    if (argc > 7) {
        max_value_iterations = atoi(argv[7]);
    }

    

    std::vector<real> rewards(horizon);
    for (int t=0; t<horizon; t++) {
        rewards[t] = 0.0;
    }


    
    real average_oracle_return = 0.0;
    
    RandomNumberFile rng("./dat/r1e7.bin");

    for (int experiment=0; experiment<n_experiments; experiment++) {
        real actual_probability = rng.uniform();//urandom(0,1); //true_random(false); //urandom(0, 1);
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

            int action = MakeDecision(expansion_method,
                                      n_states,
                                      n_actions,
                                      belief,
                                      state,
                                      gamma,
                                      n_iter,
                                      verbose,
                                      max_value_iterations,
                                      fout);


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

            if (state==1) {
                break;
            }
        }
        
        real oracle_return = 0.0;
        if (actual_probability > 0.5) {
            oracle_return = (2.0 * actual_probability - 1.0)/(1.0 - gamma);
        }
        average_oracle_return += oracle_return;
    }

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
    discounted_reward *= inv_exp;
    total_reward *= inv_exp;



    printf("%d %f %d %f %f %f\n",
           expansion_method,
           gamma,
           n_iter,
           total_reward,
           discounted_reward,
           average_oracle_return);
    return 0;
}

int MakeDecision(ExpansionMethod expansion_method,
                 int n_states,
                 int n_actions,
                 SimpleBelief prior,
                 int state,
                 real gamma,
                 int n_iter,
                 int verbose,
                 int max_value_iterations,
                 FILE* fout)
{
    BeliefTree<SimpleBelief> tree(prior, state, n_states, n_actions);
    
    for (int iter=0; iter<n_iter; iter++) {
        std::vector<BeliefTree<SimpleBelief>::Node*> node_set = tree.getNodes();
        int n_edge_nodes = 0;
        int node_index = -1;
        std::vector<int> leaf_nodes;
        for (uint i=0; i<node_set.size(); ++i) {
            if (node_set[i]->outs.size()==0) {
                leaf_nodes.push_back(i);
            } else {
                n_edge_nodes++;
            }
        }
        
        int n_leaf_nodes = (int) leaf_nodes.size();
        
        if (expansion_method == SerialExpansion) {
            node_index = leaf_nodes[0];
        } else if (expansion_method == RandomExpansion) {
            int X =  floor(urandom()*((real) leaf_nodes.size()));
            node_index = leaf_nodes[X];
        } else if (expansion_method == HighestMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = node->belief.getGreedyReturn(node->state, gamma);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == HighestDiscountedMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = ((real) node->depth) * log(gamma)
                    + log(node->belief.getGreedyReturn(node->state, gamma));
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = node->belief.sampleReturn(node->state, gamma);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = ((real) node->depth) * log(gamma)
                    + log(node->belief.sampleReturn(node->state, gamma));
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonBound) {
            std::vector<real> U(leaf_nodes.size());
            //std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
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
                BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  ((real) node->depth) * log(gamma) + log(Ui);
                //L[i] = Li;
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else {
            std::cerr << "Unknown method " << expansion_method << std::endl;
            exit(-1);
        }

        if (verbose >= 100) {
            std::cout << n_leaf_nodes << " leaf nodes, "
                      << n_edge_nodes << " edge nodes, "
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
    value_iteration.ComputeStateValues(0.00001, max_value_iterations);
    end_time = GetCPU();
    if (verbose) {
        std::cout << "# CPU " << end_time - start_time << std::endl;
    }

    if (verbose >= 60) {
        for (int s=0; s<mdp.getNStates(); s++) {
            std::cout << "V[" << s << "]"
                      << " = " << value_iteration.getValue(s)
                      << std::endl;
        }
    }

    if (fout) {
        for (int s=0; s<mdp.getNStates(); s++) {
            //BeliefTree<SimpleBelief>::Node* node = tree->nodes[s];
            real alpha = 0.0;
            real beta = 0.0;
            if (s < mdp.getNStates() - 1) {
                BetaDistribution belief = tree.nodes[s]->belief.getPrior();
                alpha = belief.alpha;
                beta = belief.beta;
            }
            fprintf (fout,
                     "s%d [label = \"%.2f %d %d\"];\n",
                     s, value_iteration.getValue(s),
                     (int) alpha, (int) beta);
            //belief.alpha,
            //       belief.beta);
        }
    }
#endif

    start_time = GetCPU();
    value_iteration.ComputeStateActionValues(0.00001, max_value_iterations);
    end_time = GetCPU();
    if (verbose) {
        std::cout << "# CPU " << end_time - start_time << std::endl;
    }

    if (verbose >= 10) {
        for (int s=0; s<mdp.getNStates(); s++) {
            for (int a=0; a<mdp.getNActions(); a++) {
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
