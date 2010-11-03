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

#ifndef COVER_TREE_H
#define COVER_TREE_H

#include "real.h"
#include "Vector.h"

#undef DEBUG_COVER_TREE
#undef DEBUG_COVER_TREE_NN

/** A cover tree.

    Implements the algorithm "Cover trees for Nearest Neighbor",
    Beygelzimer, Kakade, Langford. Based on the technical report "Fast
    Nearest Neighbors" of Thomas Kollar.

	The tree is constructed from a set \f$S\f$ of points, a metric
	\f$\rho\f$ and a constant \f$c > 1\f$. We write \f$d_{i,j} =
	\rho(x_i, x_j)\f$ for the distance between the i-th and j-th
	points. Let \f$S_n\f$ denote the set of nodes at level \f$n\f$ of
	the tree and \f$C(i)\f$ the set of children of the i-th node. The
	tree has the following properties.

	1. If \f$i \in S_n\f$ then \f$i \in S_{n-1}\f$.

	2. \f$S_\infty = \fS$.

	3. For any \f$i, j \in S_k\f$, \f$d_{i,j} > c^d\f$.
	
	4. If \f$i \in S_n\f$ and \f$j \in C(i)\f$, 
    
 */
class CoverTree
{
public:
	/// The raw metric used between points.
	static const real metric(const Vector& x, const Vector& y)
	{
		return L1Norm(&x, &y);
	}

    /// This simply is a node
    struct Node
    {
        Vector point; ///< The location of the point.
        std::vector<Node*> children; ///< Pointer to children
        int level; ///< Level in the tree
		int children_level; ///< Level of children

		/// Constructor needs a point and a level
        Node (const Vector& point_, const int level_);

		/// Destructor
        ~Node();

		/// Metric
		const real distanceTo (const Vector& x) const
		{
			return CoverTree::metric(x, point);
		}

		/// Insert a new point at the given level, as a child of this node
        const Node* Insert(const Vector& new_point, const int level);
        
        /// Find nearest neighbour of the node
        std::pair<const CoverTree::Node*, real> NearestNeighbour(const Vector& query, const real distance) const;

 		/// The number of children
        const int Size() const
        {
            return children.size();
        }
		
		/// A lower bound on the level of children
		const int ChildrenLevel() const
		{
			return children_level;
		}

		/// Display the tree in textual format
		void Show() const;

		/// Display the tree in .dot format
		void Show(FILE* fout) const;
    };

    /** Cover set.
		
		This is primarily a simple encapsulation of a vector of Node*.
		On the other hand, it also provides some set semantics:
		a node value may not be duplicated within a set.
	 */
    struct CoverSet
    {
        std::vector<Node*> nodes;
        std::vector<real> distances;
        const int Size() const
        {
            return nodes.size();
        }
        void Insert(Node* node, real distance)
        {
			int N = Size();
			for (int i=0; i<N; ++i) {
				if (nodes[i] == node) {
					return;
				}
			}
            nodes.push_back(node);
            distances.push_back(distance);
        }
        const int NChildren(const int i) const
        {
            return nodes[i]->Size();
        }
		void Show() const 
		{
			for (int i=0; i<Size(); ++i) {
				printf ("%d : ", i);
				nodes[i]->Show();
			}
		}
    };
	
	const real metric(const CoverSet& Q, const Vector& p) const;
	const Node* Insert(const Vector& new_point, const CoverSet& Q_i, const int level);
	const Node* Insert(const Vector& new_point);


	const Node* NearestNeighbour(const Vector& query_point) const;
	bool Check() const;
	void Show() const;
    CoverTree();
    ~CoverTree();
protected:
	bool Check(const CoverSet& parents, const int level) const;
	real Separation(const CoverSet& Q) const;
	int tree_level;
    Node* root;
};


#endif
