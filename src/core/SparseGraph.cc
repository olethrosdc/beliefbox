// -*- Mode: c++ -*-
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: Graph.c,v 1.4 2005/10/20 14:45:27 olethros Exp $
 
 
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/**
   \file SparseGraph.cc
   
   \brief Graph functions for sparse graph.
   
   This file provides a number of specialised operations for sparse graphs.
*/

#include "SparseGraph.h"
#include <iostream>
#include <list>
#include <exception>
#include <stdexcept>
#include <cassert>

/// Constructor, initialises the data structures.
SparseGraph::SparseGraph(int N, bool directional) : Graph(N, directional)
{
    parents.resize(N);
    children.resize(N);
}

/** Add an edge to the graph.
 * 
 * Add an edge \c e to the graph.  An old edge connecting the same two nodes
 * will be replaced if \c clear is true.  (Defaults to false).
 */
bool SparseGraph::AddEdge(Edge e, bool clear)
{
    if (!(e.src>=0 && e.src < N)) {
        std::cerr <<e.src << " " <<  e.src << " " << N << std::endl;
    }
    if (!(e.dst>=0 && e.dst < N)) {
        std::cerr <<e.dst << " " <<  e.dst << " " << N << std::endl;
    }
    assert(e.dst>=0 && e.dst < N);
    assert(e.src>=0 && e.src < N);
    assert(e.dst>=0 && e.dst < N);
    HalfEdge parent;
    HalfEdge child;
    parent.node = e.src;
    parent.w = e.w;
    child.node = e.dst;
    child.w = e.w;
    if (edge(e.src, e.dst)) {
        if (clear) {
            getParentHalfEdge(e.src, e.dst).w = e.w;
            getChildHalfEdge(e.src, e.dst).w = e.w;
        }
    } else {
        parents[e.dst].push_back(parent);
        children[e.src].push_back(child);
    }
#if GRAPH_DEBUG_LEVEL > 99
    printf ("adding edge: %d -> %d\n", e.src, e.dst);
#endif
    return true;
}

/// Recursive implementation of Dijkstra's algorithm.
///
/// Rewritten to avoid the loop.
bool SparseGraph::rCalculateDistance_iter (real* dist, int j, int* C)
{
    real current_dist = dist[j];
    HalfEdgeList& parent_list = parents[j];
    for (HalfEdgeListIterator i=parent_list.begin();
         i!=parent_list.end();
         ++i) {
        int n = i->node;
        if (C==NULL || C[n]==0) {
            if (n!=j) {
                real d = i->w + current_dist;
                if (dist[n]<0.0 || d<dist[n]) {
					dist[n] = d;
					rCalculateDistance_iter(dist, n, C);
				}
            }
        }
        
    }
	return true;
}


bool SparseGraph::edge (int src, int dst)
{
    assert(src>=0 && src < N);
    assert(dst>=0 && dst < N);
    HalfEdgeList& child_list = children[src];
    for (HalfEdgeListIterator i=child_list.begin();
         i!=child_list.end();
         ++i) {
        if (dst == i->node) {
            return true;
        }
    }
    return false;
}



real SparseGraph::distance (int src, int dst)
{
    assert(src>=0 && src < N);
    assert(dst>=0 && dst < N);
    HalfEdgeList& child_list = children[src];
    for (HalfEdgeListIterator i=child_list.begin();
         i!=child_list.end();
         ++i) {
        if (dst == i->node) {
            return i->w;
        }
    }
    return -1;
}

HalfEdge& SparseGraph::getParentHalfEdge (int src, int dst)
{
    assert(src>=0 && src < N);
    assert(dst>=0 && dst < N);
    HalfEdgeList& child_list = children[src];
    for (HalfEdgeListIterator i = child_list.begin();
         i!=child_list.end();
         ++i) {
        if (dst == i->node) {
            return *i;
        }
    }
    throw std::runtime_error("There should've been a half edge here");
}

HalfEdge& SparseGraph::getChildHalfEdge (int src, int dst)
{
    assert(src>=0 && src < N);
    assert(dst>=0 && dst < N);
    HalfEdgeList& parent_list = parents[dst];
    for (HalfEdgeListIterator i=parent_list.begin();
         i!=parent_list.end();
         ++i) {
        if (src == i->node) {
            return *i;
        }
    }
    throw std::runtime_error("There should've been a half edge here");
}

HalfEdgeListIterator SparseGraph::getFirstParent(int node)
{
    return parents[node].begin();
}

HalfEdgeListIterator SparseGraph::getFirstChild(int node)
{
    return children[node].begin();
}

int SparseGraph::n_parents(int node)
{
    return parents[node].size();
}

int SparseGraph::n_children(int node)
{
    return children[node].size();
}

/// Detecting cylces
bool SparseGraph::hasCycles_iter (std::vector<bool>& mark, int n)
{
    mark[n] = true;
    HalfEdgeList& child_list = children[n];
    for (HalfEdgeListIterator i=child_list.begin();
         i!=child_list.end();
         ++i) {
        int child = i->node;
        // for any child check to see if it's already marked
        // or (if not) whether there are subsequent cycles
        if (mark[child]) {
#if GRAPH_DEBUG_LEVEL > 50
            printf ("cycle: %d -> %d\n", n, child);
#endif
            return true; 
        } else if (hasCycles_iter(mark, child)) { 
            return true;
        }
    }
    // if no children are marked or have cycles
	return false;
}
