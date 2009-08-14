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

#ifndef KD_TREE_H
#define KD_TREE_H

#include "Vector.h"

#include <list>

class KDNode
{
 public:
    Vector box_sup;
    Vector box_inf;
    Vector c;
    int a;
    KDNode* lower;
    KDNode* upper;
    
    KDNode(Vector& c_, int a_, Vector& inf, Vector& sup) : c(c_), a(a_), lower(NULL), upper(NULL)
    {
        assert(a >= 0 && a < c.Size());
        box_inf = inf;
        box_sup = sup;
        // nothing else to init.
    }

    KDNode* AddVector(Vector& x, Vector& inf, Vector& sup);
    void NearestNeighbour(Vector& x, KDNode*& nearest, real& dist);

};

/** KD Tree

    The algorithm idea is very simple.
    
    Traverse the tree. At each node \f$i\f$ of the tree, find whether
    we are in the upper or lower half of the limits, so whether \f$x
    \in L_i\f$ or \f$x \in U_i\f$. Let us call this subset \f$X'_i\f$.
    
    If there already exists a node \f$j\f$ there then we repeat.
    Otherwise, we add a new node.
   
 */
class KDTree
{
protected:
    int n_dimensions;
    Vector box_sup;
    Vector box_inf;
    KDNode* root;
    std::vector<KDNode*> node_list;
public:	
    KDTree(int n);
    ~KDTree();
    void AddVector(Vector& x);
    void Show();
    KDNode* FindNearestNeighbourLinear(Vector& x);
    KDNode* FindNearestNeighbour(Vector& x);
};


#endif
