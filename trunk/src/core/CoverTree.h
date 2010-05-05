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
            : point(point_), level(level_)
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
			int N = Size();
			for (int i=0; i<N; ++i) {
				if (nodes[i] == node) {
					//printf("Weird!\n");
					return;
				}
			}
            nodes.push_back(node);
        }
        int NChildren(int i)
        {
            return nodes[i]->Size();
        }
    };

	real metric(CoverSet& Q, Vector& p)
	{
		real D = INF;
		for (int i=0; i<Q.Size(); ++i) {
			real d_i = metric(Q.nodes[i]->point, p);
			if (d_i < D) {
				D = d_i;
			}
		}
		return D;
	}

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
                real dist_i = metric(new_point, node->point);
                
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
        bool found = Insert(new_point, Q_next, level - 1);
        real distance = INF;
        for (int k=0; k<Q_i.Size(); ++k) {
            Node* node = Q_i.nodes[k];
            real dist_k = metric(new_point, node->point);
            if (dist_k < distance) {
                distance = dist_k; 
                closest_node = node;
                if (distance <= separation) {
                    break;
                }
            }
        }

        if (found && distance <= separation) {
            closest_node->Insert(new_point, level - 1);
#ifdef DEBUG_COVER_TREE
            printf("Inserted at level %d\n", level - 1);
#endif
            return false;
        }
        return found;
   }


	/// Insert a new point in the tree
    void Insert(Vector& new_point)
    {
        if (!root) {
#ifdef DEBUG_COVER_TREE
            printf("Adding root at:");
            new_point.print(stdout);
            printf("\n");
#endif
            root = new Node(new_point, INF);
            return;
        }
#ifdef DEBUG_COVER_TREE
		printf("Trying to add new point\n");
#endif
		real distance = metric(new_point, root->point);
		int level = (int) ceil(log(distance) / log(2));
        CoverSet Q;
        Q.Insert(root);
        Insert(new_point, Q, level);
    }


    /// Find the nearest node
    Node* NearestNeighbour(Vector& query_point, CoverSet& Q_i, int level)
    {
		CoverSet Q_next;
		CoverSet Q;
        // go through all the children and only add them if they are close
		int next_level = level - 1;

		// Get Q
        for (int k=0; k<Q_i.Size(); ++k) {
            int n_children = Q_i.NChildren(k);
#ifdef DEBUG_COVER_TREE
			printf("Node: %d with %d children:", k, n_children); Q_i.nodes[k]->point.print(stdout);
#endif
            for (int j= -1; j<n_children; ++j) {
                Node* node;
                if (j >= 0) {
                    node = Q_i.nodes[k]->children[j];
                } else {
                    node = Q_i.nodes[k];
                }
				Q.Insert(node);
#ifdef DEBUG_COVER_TREE
				printf("Q: [l:%d] ", node->level); node->point.print(stdout);
#endif
			}
		}

		real dist_Q_p = metric(Q, query_point);
		real separation = dist_Q_p + pow(2, next_level);
		Node* found_node = NULL;
		real min_dist = separation * 2;
		bool close_children = false;
		for (int i=0; i<Q.Size(); ++i) {
			Node* node = Q.nodes[i];
			real dist_i = metric(query_point, node->point);
			if (dist_i <= separation) {
				Q_next.Insert(node);
				close_children = true;
			} 
			if (dist_i < min_dist) {
				min_dist = dist_i;
				found_node = node;
			}
		}
		if (!close_children) {
#ifdef DEBUG_COVER_TREE
            printf("Found node at level %d\n", level);
#endif
			return found_node;
		} else {
#ifdef DEBUG_COVER_TREE
			printf("Distance: %f, jumping down to level %d\n", min_dist, next_level);
#endif
			return NearestNeighbour(query_point, Q_next, next_level);
		}
   }

	/// FInd the nearest neighbour in the tree
    Node* NearestNeighbour(Vector& query_point)
    {
        if (!root) {
            return NULL;
        }
#ifdef DEBUG_COVER_TREE
        printf("Query: "); query_point.print(stdout);
#endif
		real distance = metric(query_point, root->point);
		int level = (int) ceil(log(distance) / log(2));
        CoverSet Q;
        Q.Insert(root);
        return NearestNeighbour(query_point, Q, level);
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
