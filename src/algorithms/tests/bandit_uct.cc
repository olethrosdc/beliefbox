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


//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


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

    RandomNumberFile rng("./dat/r1e7.bin");
    int n_states  = 0;

    // perform experiments
    for (int experiment=0; experiment<n_experiments; experiment++) {
        
        std::vector<real> Er(n_actions);

        for (int i=0; i<n_actions; i++) {
            Er[i] = rng.uniform();
        }
        
        // initial state and belief
        int state = 0;
        BanditBelief belief(n_actions, alpha, beta);
                                                
        // looop over time
        for (int t=0; t<horizon; t++) {
            // write graph to a file if we are only doing one experiment
            FILE* fout = NULL;
            if (n_experiments == 1) {
                char buffer[1024];
                sprintf(buffer, "test%d.dot", t);
                fout = fopen (buffer, "w");
                fprintf (fout, "digraph Lookahead {\n");
                fprintf (fout, "ranksep=2; rankdir=LR; \n");
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

    // show statistics
    printf("%d %d %f %d %f %f %f %f\n",
           method,
           n_iter,
           gamma,
           n_actions,
           total_reward,
           discounted_reward,
           average_oracle_reward,
           average_oracle_return);
    return 0;
}

#if 0
static bool sort_using_greater_than(real u, real v)
{
    return u > v;
}
#endif

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
    BeliefTree<BanditBelief> tree(prior, state, n_states, n_actions);
    std::vector<BeliefTree<BanditBelief>::Node*> node_set = tree.getNodes();

    for (int iter=0; iter<n_iter; iter++) {
        std::vector<BeliefTree<BanditBelief>::Node*> node_set = tree.getNodes();
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
            int X =  (int) floor(urandom()*((real) leaf_nodes.size()));
            node_index = leaf_nodes[X];
        } else if (expansion_method == HighestMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<BanditBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = node->belief.getGreedyReturn(node->state, gamma);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == HighestDiscountedMeanValue) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<BanditBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = ((real) node->depth) * log(gamma)
                    + log(node->belief.getGreedyReturn(node->state, gamma));
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<BanditBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = node->belief.sampleReturn(node->state, gamma);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedThompsonSampling) {
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<BanditBelief>::Node* node = node_set[leaf_nodes[i]];
                U[i] = ((real) node->depth) * log(gamma)
                    + log(node->belief.sampleReturn(node->state, gamma));
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == ThompsonBound) {
            std::vector<real> U(leaf_nodes.size());
            //std::vector<real> L(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                BeliefTree<BanditBelief>::Node* node = node_set[leaf_nodes[i]];
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
                BeliefTree<BanditBelief>::Node* node = node_set[leaf_nodes[i]];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  ((real) node->depth) * log(gamma) + log(Ui);
                //L[i] = Li;
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == HighProbabilityBoundOnly) {
            // returns the maximum of the current upper bounds
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTree<BanditBelief>::Node* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                U[i] =  log(Ui);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedHighProbabilityBoundOnly) {
            // returns the maximum of the current upper bounds
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTree<BanditBelief>::Node* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                U[i] =  ((real) node->depth) * log(gamma) + log(Ui);
            }
            node_index = leaf_nodes[ArgMax(U)]; 
        } else if (expansion_method == HighProbabilityBound) {
            // returns the maximum of the current upper bounds or the lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTree<BanditBelief>::Node* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  log(Ui);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedHighProbabilityBound) {
            // returns the maximum of the current upper bounds or the lower bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTree<BanditBelief>::Node* node = node_set[n];
                node_set[n]->U.push_back(node->belief.sampleReturn(node->state, gamma));
                real Ui = Max(node_set[n]->U);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  ((real) node->depth) * log(gamma) + log(Ui);
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == MeanHighProbabilityBound) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTree<BanditBelief>::Node* node = node_set[n];
                std::vector<real> &Ub = node_set[n]->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Li > Ui) {
                    Ui = Li;
                }
                U[i] = Ui;
            }
            node_index = leaf_nodes[ArgMax(U)];
        } else if (expansion_method == DiscountedMeanHighProbabilityBound) {
            // returns the mean high probability bound
            std::vector<real> U(leaf_nodes.size());
            for (int i=0; i<n_leaf_nodes; i++) {
                int n = leaf_nodes[i];
                BeliefTree<BanditBelief>::Node* node = node_set[n];
                std::vector<real> &Ub = node_set[n]->U;
                int n_samples=1;
                for (int k=0; k<n_samples; k++) {
                    Ub.push_back(node->belief.sampleReturn(node->state, gamma));
                }
                real Ui = Mean(Ub);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Li > Ui) {
                    Ui = Li;
                }
                U[i] =  ((real) node->depth) * log(gamma) + log(Ui);
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


    if (fout) {
        value_iteration.ComputeStateValues(0.00001, max_value_iterations);
        if (verbose >= 60) {
            for (int s=0; s<mdp.GetNStates(); s++) {
                std::cout << "V[" << s << "]"
                          << " = " << value_iteration.getValue(s)
                          << std::endl;
            }
        }

        for (int s=0; s<mdp.GetNStates(); s++) {
            //BeliefTree<BanditBelief>::Node* node = tree->nodes[s];
            //real alpha = 0.0;
            //real beta = 0.0;
            //            if (s < mdp.GetNStates() - 1) {
            //  std::vector<BetaDistribution*> belief = tree.nodes[s]->belief.getPrior();
            //   alpha = belief.alpha;
            //   beta = belief.beta;
            //}
            fprintf (fout,
                     "s%d [label = \"%.2f\"];\n",
                     s, value_iteration.getValue(s));
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
