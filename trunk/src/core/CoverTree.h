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

#define DEBUG_COVER_TREE
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
        int level;
        Node (Vector& point_, int level_)
            : point(point_, int level_)
        {
        }
        ~Node()
        {
            for (uint i=0; i<children.size(); ++i) {
                delete children[i];
            }
        }
        void Insert(Vector& new_point, int level)
        {
#ifdef DEBUG_COVER_TREE
            printf("New child for: ");
            point.print(stdout);
            printf("at: ");
            new_point.print(stdout);
#endif            
            Node* node = new Node(new_point, level);
            children.push_back(node);
        }
        int Size()
        {
            return children.size();
        }
    };

    /// Cover set
    struct CoverSet
    {
        std::vector<Node*> nodes;
        int Size()
        {
            return nodes.size();
        }
        void Insert(Node* node)
        {
            nodes.push_back(node);
        }
        int NChildren(int i)
        {
            return nodes[i]->Size();
        }
    };

    /// Insert
    bool Insert(Vector& new_point, CoverSet& Q_i, int level)
    {
        Node* closest_node = NULL;
		// Check if d(p, Q) > 2^level
		real log_separation = level * log(2);
		real separation = exp(log_separation);
		bool separated = true;
		CoverSet Q_next = Q_i;

        // go through all the children and only add them if they are close
        for (int k=0; k<Q_i.Size(); ++k) {
            int n_children = Q_i.NChildren(k);
            for (int j= -1; j<n_children; ++j) {
                Node* node;
                if (j >= 0) {
                    node = Q_i.nodes[k]->children[j];
                    if (node->level != level - 1) {
                        continue;
                    }
                } else {
                    node = Q_i.nodes[k];
                }
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
        real distance = INF;
        for (int k=0; k<Q_i.Size(); ++k) {
            Node* node = Q_i.nodes[k];
            real dist_k = metric(p, node->point);
            if (dist_k < distance) {
                distance = dist_k; 
                closest_node = node;
                if (distance <= separation) {
                    break;
                }
            }
        }

        if (found && distance <= separation) {
            closest_node->Insert(p);
#ifdef DEBUG_COVER_TREE
            printf("Inserted at level %d\n", level - 1);
#endif
            return false;
        }
        return found;
   }

    void Insert(Vector& p)
    {
        if (!root) {
#ifdef DEBUG_COVER_TREE
            printf("Adding root at:");
            p.print(stdout);
            printf("\n");
#endif
            root = new Node(p);
            return;
        }
#ifdef DEBUG_COVER_TREE
            printf("Trying to add new point\n");
#endif
            //real distance = metric(p, root->point);
        int level = 0;//(int) ceil(log(distance) / log(2));
        CoverSet Q;
        Q.Insert(root);
        Insert(p, Q, level);
    }

    void Show()
    {
        
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
