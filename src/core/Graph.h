/* -*- Mode: C++; -*- */
/* VER: $Id: StateActionPolicy.h,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GRAPH_H
#define GRAPH_H

#include <cstdlib>
#include <list>
#include <vector>
#include <set>
#include "real.h"

/**
   \defgroup GraphTheory Graph theory

   There are three possible methods to make graphs.

   1. A connectivity matrix.  This is very simple but takes up N^2
   space.  This is implemented as the class ConnectivityMatrix.
   Finding all neighoubrs for a node is O(N).  Finding whether i,j
   are connected is O(1).

   2. A list of edges for each node.  This is somewhat more
   complicated but takes up O(N) space.  We can use a single list
   and just have a pointer for every node.  This is implemented as
   the class SparseGraph.  Finding whether i,j are connected is
   O(N).  Insertion is O(1).  Finding all neighbours for a node is
   O(N). 
*/
/*@{*/


/** Structure holding an edge
 */
struct Edge
{
    int src;
    int dst;
    real w;
    Edge (int src_, int dst_, real w_ = 1.0)
    : src(src_), dst(dst_), w(w_)
    {
    }
};

/** Structure holding a half-edge
 */
struct HalfEdge
{
    int node;
    real w;
};


#define GRAPH_DEBUG_LEVEL 50



/// A general graph class.  This class provides an interface for
/// classes implementing graph semantics. Most of the non-pure virtual
/// functions do not need to be over-ridden, unless a sub-class has a
/// special fast implementation for one of them.
class Graph
{
protected:
    int N; ///< number of nodes
    bool directional; ///< whether it is directional
    virtual bool rCalculateDistance_iter(real* dist, int j, int* C=NULL);
    virtual bool hasCycles_iter (std::vector<bool>& mark, int n);
public:
    Graph(int N, bool directional);
    virtual ~Graph();
    virtual bool edge (int src, int dst) = 0; ///< Returns true if there is an edge from \c src to \c dst.
    virtual real distance (int src, int dst) = 0; ///< Distance between two nodes. Note that in some applications this would have the meaning of an edge weight rather than a distance.
    int n_nodes() {return N;}
    bool is_directional() {return directional;}
    virtual int n_out_edges (int i);
    virtual int n_in_edges (int i);
    virtual int n_edges (int i);
    virtual int kConnectivity(); 
    virtual int kConnectivity(int i, int* C = NULL); 
    virtual int kConnectivity(int i, int j, int* C = NULL); 
    virtual bool CalculateDistances (real* dist, int j, int* C=NULL);
    virtual bool isReachable (int j, int* C=NULL);
    virtual bool MarkShortestPath (real* dist, int* C, int* P, int i, int j);
    virtual bool hasCycles ();
};
/*@}*/




#endif
