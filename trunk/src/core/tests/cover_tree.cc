/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "EasyClock.h"
#include "CoverTree.h"
#include "KDTree.h"
#include "Random.h"
#include "Vector.h"
#include "NormalDistribution.h"
#include <vector>
#include <memory>

void test_cover_tree_insertion(CoverTree& tree, std::vector<Vector>& X)
{
    int n_points = X.size();

    for (int i=0; i<n_points; i++) {
        //printf ("Adding : "); X[i].print(stdout);
        CoverTree::Node* new_node = tree.Insert(X[i]);
        assert(new_node);
    }    
}

int test_cover_tree_query(CoverTree& tree, std::vector<Vector>& Q)
{
    int n_points = Q.size();
    int n_errors = 0;
    for (int i=0; i<n_points; ++i) {
        const CoverTree::Node* node = tree.NearestNeighbour(Q[i]);
        if (!node) {
            n_errors++;
            printf ("Null node!\n");
        }
    }
    return n_errors;
}

void test_kd_tree_insertion(KDTree<void>& tree, std::vector<Vector>& X)
{
    int n_points = X.size();

    for (int i=0; i<n_points; i++) {
        //printf ("Adding : "); X[i].print(stdout);
        tree.AddVectorObject(X[i], NULL);
    }    
}

int test_kd_tree_query(KDTree<void>& tree, std::vector<Vector>& Q)
{
    int n_points = Q.size();
    int n_errors = 0;
    for (int i=0; i<n_points; ++i) {
        KDNode* node = tree.FindNearestNeighbour(Q[i]);
        if (!node) {
            n_errors++;
        }
    }
    return n_errors;
}

int check_cover_tree_query(CoverTree& tree,
                           std::vector<Vector>& X,
                           std::vector<Vector>& Q)
{
    int n_errors = 0;
    int n_stored_points = X.size();
    int n_test_points = Q.size();
    for (int i=0; i<n_test_points; ++i) {
        const CoverTree::Node* node = tree.NearestNeighbour(Q[i]);

        if (node) {
            Vector best_point = node->point;
            real dist = INF;
            int arg_min = -1;
            for (int k=0; k<n_stored_points; ++k) {
                real dist_k = tree.metric(X[k], Q[i]);
                if (dist_k < dist) {
                    dist = dist_k;
                    arg_min = k;
                }
            }
        
            real error = tree.metric(best_point, X[arg_min]);
        
            if (error > 0.000001) {
                printf("## Distance is too big: %f ##########\n", error);
                printf("Query: "); Q[i].print(stdout);
                printf("Found: (%f) ", tree.metric(Q[i], node->point)); best_point.print(stdout);
                printf("Best: (%f) ", dist); X[arg_min].print(stdout);
                Vector Q_alt = X[arg_min];
                const CoverTree::Node* node_alt = tree.NearestNeighbour(Q_alt);
                printf("Alternative result: "); node_alt->point.print(stdout);
                n_errors++;
            }
        } else {
            printf ("NULL\n");
            n_errors++;
        }
    }
    return n_errors;
}


int main(int argc, char** argv)
{

	if (argc != 5) {
		fprintf(stderr, "Usage: cover_tree points dimensions n_trials c\n");
		exit(-1);
	}
    int n_points = atoi(argv[1]);
    int n_dimensions = atoi(argv[2]);
    int n_trials = atoi(argv[3]);
    real c = atof(argv[4]);
    printf ("#  %d dataset points, %d test points, %d dimensions\n",
            n_points, n_trials, n_dimensions);
    printf("# Generating data\n");

    std::auto_ptr<Distribution> distribution(new NormalDistribution);
    Vector Transform(n_dimensions);
    for (int j=0; j<n_dimensions; j++) {
        Transform(j) = distribution->generate();
        if (j) {
            Transform(j) = 0;
        }
    }
    std::vector<Vector> X(n_points);
    for (int i=0; i<n_points; i++) {
        X[i].Resize(n_dimensions);
        for (int j=0; j<n_dimensions; j++) {
            real Z = distribution->generate();
            X[i][j] = Z * Transform(j);
        }
    }

    std::vector<Vector> Q(n_trials);
    for (int i=0; i<n_trials; i++) {
        Q[i].Resize(n_dimensions);
        for (int j=0; j<n_dimensions; j++) {
            real Z = distribution->generate();
            Q[i][j] = Z * Transform(j);
        }
    }
    
    
    printf ("# insertion_time, query time, total time\n");
    if (1) {
        printf ("# Testing cover tree\n");
        CoverTree cover_tree(c);
        double timer_start = GetCPU();
        test_cover_tree_insertion(cover_tree, X);
        double timer_mid = GetCPU();
			//cover_tree.Show();
        int n_cover_tree_failures = test_cover_tree_query(cover_tree, Q);
        double timer_end = GetCPU();
        printf ("%f %f %f # COVER TREE\n",
                timer_mid - timer_start,
                timer_end - timer_mid,
                timer_end - timer_start);
        //n_cover_tree_failures += check_cover_tree_query(cover_tree, X, Q);
        printf ("Errors: %d\n", n_cover_tree_failures);
    }

    
    if (1) {
        printf ("# Testing KD tree\n");
        KDTree<void> kd_tree(n_dimensions);
        
        double timer_start = GetCPU();

        test_kd_tree_insertion(kd_tree, X);

        double timer_mid = GetCPU();
        int n_kd_tree_failures = test_kd_tree_query(kd_tree, Q);
        
        double timer_end = GetCPU();
        
        printf ("%f %f %f # KD TREE\n",
                timer_mid - timer_start,
                timer_end - timer_mid,
                timer_end - timer_start);

    }



    return 0;
}

#endif
