/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "KDTree.h"
#include "Random.h"
#include "Vector.h"
#include <vector>

bool kd_tree_test(int n_points, int n_dimensions) 
{
    printf ("Testing with %d points and %d dimensions\n", n_points, n_dimensions);
    std::vector<Vector> X(n_points);

    for (int i=0; i<n_points; i++) {
        X[i].Resize(n_dimensions);
        for (int j=0; j<n_dimensions; j++) {
            X[i][j] = urandom();
        }
    }
    
    KDTree<int> tree(n_dimensions);
    std::vector<int> number(n_points);
    for (int i=0; i<n_points; i++) {
        number[i] = i;
        tree.AddVectorObject(X[i], &number[i]);
    }    

    //tree.Show();
    // First find 1 nearest neigbour
    for (int i=0; i<n_points; i++) {
        Vector Z = X[i];
        KDNode* node = tree.FindNearestNeighbourLinear(Z);
        KDNode* node2 = tree.FindNearestNeighbour(Z);
        OrderedFixedList<KDNode> knn_list = tree.FindKNearestNeighboursLinear(Z, 1);

        KDNode* node3 = knn_list.S.front().second;
        OrderedFixedList<KDNode> knn_list2 = tree.FindKNearestNeighbours(Z, 1);

        KDNode* node4 = knn_list2.S.front().second;
        if (node != node2 || node != node3 || node != node4) {
            printf ("MISMATCH ");
            printf ("dist: %f %f\n", L1Norm(&Z, &node->c), L1Norm(&Z, &node2->c));
            return false;
        }
    }        

    for (int i=0; i<n_points; i++) {
        Vector Z(n_dimensions);
        for (int j=0; j<n_dimensions; ++j) {
            Z[j] = urandom();
        }
        KDNode* node = tree.FindNearestNeighbourLinear(Z);
        KDNode* node2 = tree.FindNearestNeighbour(Z);
        int K = (int) ceil(urandom(1,10));
        OrderedFixedList<KDNode> knn_list = tree.FindKNearestNeighboursLinear(Z, K);
        OrderedFixedList<KDNode> knn_list2 = tree.FindKNearestNeighbours(Z, K);
        
        KDNode* node3 = knn_list.S.front().second;
        KDNode* node4 = knn_list2.S.front().second;
        if (node != node2 || node != node3 || node != node4) {
            printf ("MISMATCH ");
            printf ("dist: %f %f\n", L1Norm(&Z, &node->c), L1Norm(&Z, &node2->c));
            return false;
        }
        
        std::list<std::pair<real, KDNode*> >::iterator s1 = knn_list.S.begin();
        std::list<std::pair<real, KDNode*> >::iterator s2 = knn_list.S.begin();
        for (int k=1; k<K; ++k, ++s1, ++s2) {
            node3 = s1->second;
            node4 = s2->second;
            if (node3 != node4) {
                printf ("MISMATCH : %f %f",
                        L1Norm(&Z, &node3->c), L1Norm(&Z, &node4->c));
                return false;
            }
        }
    }        

    printf ("%d/%d\n", tree.getNumberOfLeaves(), tree.getNumberOfNodes());

    return true;
    
}
int main(void) 
{
    int n_errors = 0;
    int n_tests = 10;
    for (int i=0; i<n_tests; ++i) {
        int n_dim = (int) ceil(urandom(1,10));
        int n_points = (int) ceil(urandom(1,1000));
        bool success = kd_tree_test(n_points, n_dim);
        if (!success) {
            n_errors ++;
        }
    }
    
    if (n_errors) {
        fprintf (stderr, "%d / %d tests failed\n", n_errors, n_tests);
        return n_errors;
    }
    printf ("Test successful\n");
    return 0;
}


#endif
