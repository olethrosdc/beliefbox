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
   \file AbstractGraph.cc
   
   \brief AbstractGraph theory functions.
   
   This file provides a number of basic operations on graphs. Some of
   them are useful for sanity checking, others for implementations of
   various algorithms.
*/

#include "AbstractGraph.h"
#include "SmartAssert.h"

AbstractGraph::AbstractGraph(int N, bool directional)
    {
        //DISABLED_ASSERT (N>0)(N);
        this->N = N;
        this->directional = directional;
#if GRAPH_DEBUG_LEVEL > 90
        if (directional) {
            printf("Making directional graph\n");
        } else {
            printf("Making undirected graph\n");
        }
#endif
    }
AbstractGraph::~AbstractGraph()
    {
#if GRAPH_DEBUG_LEVEL > 90
        printf("Destroying AbstractGraph\n");
#endif
    }

bool AbstractGraph::hasCycles_iter (std::vector<bool>& mark, int n)
{
    mark[n] = true;
    for (int i=0; i<N; ++i) {
        // for any child check to see if it's already marked
        // or (if not) whether there are subsequent cycles
        if (edge(n,i)) {
            if (mark[i]) {
                return true; 
            } else if (hasCycles_iter(mark, i)) { 
                return true;
            }
        }
    }
    // if no children are marked or have cycles
    return false;
}


/// Detecting cycles.  Go through all nodes and see if you 
/// see a cycle starting from any one node.
bool AbstractGraph::hasCycles ()
{
    std::vector<bool> mark(N);
    std::vector<bool> mark_init(N);
    for (int i=0; i<N; i++) {
        mark_init[i] = false;
    }
    for (int n=0; n<N; n++) {
        // only start from nodes which have not already been searched.
        if (!mark_init[n]) {
            for (int i=0; i<N; i++) {
                mark[i] = false;
            }
            if (hasCycles_iter(mark, n)) {
                return true;
            }
            for (int i=0; i<N; i++) {
                mark_init[i] = mark_init[i] || mark[i];
            }
        }
    }
    return false;
}


/// Number of outgoing edges from node i
int AbstractGraph::n_out_edges(int i)
{
    int cnt = 0;
    for (int j=0; j<N; j++){
        if (edge(i,j)) cnt++;
    }
    return cnt;
}

/// Number of incoming edges to node i
int AbstractGraph::n_in_edges(int i)
{
    int cnt = 0;
    for (int j=0; j<N; j++){
        if (edge(j,i)) cnt++;
    }
    return cnt;
}

/// Total number of edges for node i
int AbstractGraph::n_edges(int i)
{
    int cnt = 0;
    for (int j=0; j<N; j++){
        if (edge(j,i)) cnt++;
		if (i!=j && edge(i,j)) cnt++;
    }
    return cnt;


}


