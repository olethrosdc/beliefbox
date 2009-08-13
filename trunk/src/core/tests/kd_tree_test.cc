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

    tree.Show();

    for (int i=0; i<n_points; i++) {
        Vector Z = X[i];
        KDNode* node = tree.FindNearestNeighbourLinear(Z);
        KDNode* node2 = tree.FindNearestNeighbour(Z);
        if (node != node2) {
            printf ("MISMATCH ");
        }
        printf ("dist: %f %f\n", L1Norm(&Z, &node->c), L1Norm(&Z, &node2->c));
    }        

    return true;
    
}
int main(void) 
{
    bool success = kd_tree_test(10, 1);
    if (!success) {
        printf ("Test failed\n");
        return -1;
    }
    
    printf ("Test successful\n");
    return 0;
}


#endif
