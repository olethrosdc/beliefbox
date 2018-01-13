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
#include "OrderedFixedList.h"

#include <list>

/// Implementation of a KD tree node
class KDNode
{
 public:
    Vector box_sup; ///< upper bound 
    Vector box_inf; ///< lower bound
    Vector c; ///< center
    int a; ///< split dimension
    KDNode* lower; ///< lower child
    KDNode* upper; ///< upper child
    const void* object; ///< easiest way to associate an object
	
	/// Make a node
    KDNode(const Vector& c_, int a_, Vector& inf, Vector& sup, const void* object_) : c(c_), a(a_), lower(NULL), upper(NULL), object(object_)
    {
        assert(a >= 0 && a < c.Size());
        box_inf = inf;
        box_sup = sup;
        // nothing else to init.
    }

    KDNode* AddVector(const Vector& x, Vector& inf,  Vector& sup, const void* object);
    void NearestNeighbour(const Vector& x, KDNode*& nearest, real& dist);
    void KNearestNeighbours(const Vector& x,OrderedFixedList<KDNode>& knn_list, real& dist);
};

/** void_KD Tree

    The algorithm idea is very simple.
    
    Traverse the tree. At each node \f$i\f$ of the tree, find whether
    we are in the upper or lower half of the limits, so whether \f$x
    \in L_i\f$ or \f$x \in U_i\f$. Let us call this subset \f$X'_i\f$.
    
    If there already exists a node \f$j\f$ there then we repeat.
    Otherwise, we add a new node.

	The upper and lower half are defined slightly differently from the
	standard KD-tree.  Instead of splitting along the longest
	dimension at \f$x\f$, we split along the longest dimension at the
	centroid \f$c\f$ of the box.
   
 */
class void_KDTree
{
protected:
    int n_dimensions; ///< dimensionality of space
    Vector box_sup; ///< global upper bound
    Vector box_inf; ///< global lower bound
    KDNode* root; ///< root node
    std::vector<KDNode*> node_list; ///< contains a list of all nodes
public:	
    void_KDTree(int n); 
    virtual ~void_KDTree();
    void AddVector(const Vector& x, const void* object);
    void Show() const;
    KDNode* FindNearestNeighbourLinear(const Vector& x);
    KDNode* FindNearestNeighbour(const Vector& x);
    OrderedFixedList<KDNode> FindKNearestNeighboursLinear(const Vector& x, const int K) const;
    OrderedFixedList<KDNode> FindKNearestNeighbours(const Vector& x, const int K) const; 
    typedef std::list<std::pair<real, KDNode*> >::iterator iterator;
	/// Get number of nodes
    int getNumberOfNodes() const
    {
        return node_list.size();
    }
	/// Get number of leaves
    int getNumberOfLeaves() const
    {
        int n_leaves = 0;
        int N = node_list.size();
        for (int i=0; i<N; ++i) {
            if (node_list[i]->lower == NULL) {
                n_leaves++;
            }
            if (node_list[i]->upper == NULL) {
                n_leaves++;
            }
        }
        return n_leaves;
    }
};

/// This template makes the void* type safe.
template <typename T>
class KDTree : public void_KDTree
{
 public:
	/// Make a KD-Tree
    KDTree(int n) : void_KDTree(n)
    {
    }
	/// Find the nearest object, in linear time
    T* FindNearestObjectLinear(const Vector& x)
    {
        KDNode* node = void_KDTree::FindNearestNeighbourLinear(x);
        return (T*) node->object;
    }
	/// Find the nearest object, in logarithmic time.
    T* FindNearestObject(const Vector& x)
    {
        KDNode* node = void_KDTree::FindNearestNeighbour(x);
        return (T*) node->object;
    }
	/// Find the nearest neighbour in linear time.
    KDNode* FindNearestNeighbourLinear(const Vector& x)
    {
        return void_KDTree::FindNearestNeighbourLinear(x);

    }
	/// Find the nearest neighbour in logarithmic time.
    KDNode* FindNearestNeighbour(const Vector& x)
    {
        return void_KDTree::FindNearestNeighbour(x);
    }

    /// Find the K nearest objects in linear time
    T* FindKNearestObjectsLinear(const Vector& x, const int K)
    {
        KDNode* node = void_KDTree::FindKNearestNeighboursLinear(x, K);
        return (T*) node->object;
    }
    /// Find the K nearest objects
    T* FindKNearestObjects(const Vector& x, const int K)
    {
        KDNode* node = void_KDTree::FindKNearestNeighbours(x, K);
        return (T*) node->object;
    }
	    /// Find the K nearest objects in linear time
    OrderedFixedList<KDNode> FindKNearestNeighboursLinear(const Vector& x, const int K)
    {
        return void_KDTree::FindKNearestNeighboursLinear(x, K);

    }
    OrderedFixedList<KDNode> FindKNearestNeighbours(const Vector& x, const int K)
    {
        return void_KDTree::FindKNearestNeighbours(x, K);
    }


    T* getObject(KDNode* node)
    {
        return (T*) node->object;
    }
    
    void AddVectorObject(const Vector& x, T* object)
    {
        AddVector(x, (void*) object);
    }

    
};

#endif
