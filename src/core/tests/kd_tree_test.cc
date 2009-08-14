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
    
    KDTree tree(n_dimensions);
    for (int i=0; i<n_points; i++) {
        tree.AddVector(X[i]);
    }    

    //tree.Show();

    for (int i=0; i<n_points; i++) {
        Vector Z = X[i];
        KDNode* node = tree.FindNearestNeighbourLinear(Z);
        KDNode* node2 = tree.FindNearestNeighbour(Z);
        if (node != node2) {
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
        if (node != node2) {
            printf ("MISMATCH ");
            printf ("dist: %f %f\n", L1Norm(&Z, &node->c), L1Norm(&Z, &node2->c));
            return false;
        }
    }        

    return true;
    
}
int main(void) 
{
    int n_errors = 0;
    int n_tests = 10;
    for (int i=0; i<n_tests; ++i) {
        int n_dim = ceil(urandom(1,10));
        int n_points = ceil(urandom(1,10000));
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
