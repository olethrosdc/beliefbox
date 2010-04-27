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

/** A cover tree.

    Implements the algorithm "Cover trees for Nearest Neighbor",
    Beygelzimer, Kakade, Langford. Based on the technical report "Fast
    Nearest Neighbors" of Thomas Kollar.
    
 */
class CoverTree
{
	real metric(Vector& x, Vector& y)
	{
		return EuclideanNorm(&x, &y);
	}
    std::vector<Vector> data;
    struct Node
    {
        int level;
        real log_distance;
        real distance;
        Vector point;
        std::vector<Node*> children;
        Node (int level, Vector& point_, real log_c)
            : level(level_),
              point(point_)
        {
            log_distance = level * log_c;
            distance = exp(log_distance);
        }
        bool Contains(Vector& query) 
        {
            if (metric(query, point) <= distance) {
                return true;
            }
            return false;
        }
    };

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

    bool Insert(Vector p, CoverSet Q_i, int level)
    {
        Node* closest_node = NULL;
        real distance = INF;

		// Check if d(p, Q) > 2^level
		real log_separation = level * log(2);
		real separation = exp(log_separation);
		bool separated = false;
		CoverSet Q_next;
        for (int i=0; i<Q_i.Size(); ++i) {
            int n_children = Q_i.NChildren(i);
            for (int j=0; j<n_children; ++j) {
                Node& node = Q_i.points[i].children[j];
                real dist_i = metric(p, node->point);
                if (dist_i < distance ) {
                    distance = dist_i;
					closest_node = &node;
				}
				if (dist_i <= separation) {
					separated = true; 
					Q_next.Insert(node);
				}
            }
        } 
   }
    
};


#endif
