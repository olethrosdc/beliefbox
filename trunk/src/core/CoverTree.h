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

/** A cover tree.

    Implements the algorithm "Cover trees for Nearest Neighbor",
    Beygelzimer, Kakade, Langford. Based on the technical report "Fast
    Nearest Neighbors" of Thomas Kollar.
    
 */
class CoverTree
{
public:
	real metric(Vector& x, Vector& y)
	{
		return EuclideanNorm(&x, &y);
	}

    /// This simply is a node
    struct Node
    {
        Vector point;
        std::vector<Node*> children;
        Node (Vector& point_)
            : point(point_)
        {
        }
        ~Node()
        {
            for (uint i=0; i<children.size(); ++i) {
                delete children[i];
            }
        }
        void Insert(Vector& p)
        {
            Node* node = new Node(level - 1, p, log(2));
            children.push_back(node);
        }
    };

    /// Cover set
    struct CoverSet
    {
        std::vector<Node&> points;
        int Size()
        {
            return points.size();
        }
        void Insert(Node& point)
        {
            points.push_back(point);
        }
        int NChildren(int i)
        {
            return points[i].size();
        }
    };

    /// Insert
    bool Insert(Vector& p, CoverSet& Q_i, int level)
    {
        Node* closest_node = NULL;
        real distance = INF;
		// Check if d(p, Q) > 2^level
		real log_separation = level * log(2);
		real separation = exp(log_separation);
		bool separated = true;
		CoverSet Q_next;

        // go through all the children and only add them if they are close
        for (int i=0; i<Q_i.Size(); ++i) {
            real dist = metric(p, Q_[i].points[i].point);
            if (dist < distance_Q_i) {
                distance = dist;
                if (dist < distance ) {
                    distance = dist;
					closest_node = &node;
				}
            }
            int n_children = Q_i.NChildren(i);
            for (int j=0; j<n_children; ++j) {
                Node& node = Q_i.points[i].children[j];
                real dist_i = metric(p, node->point);
				if (dist_i <= separation) {
					separated = false; 
					Q_next.Insert(node);
				}
            }
        }
        // if all of the children are furhter away than 2^level, then
        // parent is found
        if (separated) {
            return true;
        }
        // Try and insert a new point
        bool found = Insert(p, Q_next, level - 1);
        if (found && distance <= separation) {
            closest_node->Insert(p);
            return false;
        }
        return found;
   }

    void Insert(Vector& p)
    {
        if (!root) {
            root = new Node(p);
            return;
        }
        real distance = metric(p, root->point);
        int level = (int) ceil(log(distance) / log(2));
    }

    CoverTree()
    {
        root = NULL;
    }

    ~CoverTree()
    {
        delete root;
    }
    Node* root;
};


#endif
