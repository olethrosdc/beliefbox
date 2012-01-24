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
   \file Graph.cc
   
   \brief Graph theory functions.
   
   This file provides a number of basic operations on graphs. Some of
   them are useful for sanity checking, others for implementations of
   various algorithms.
*/

#include "Graph.h"
#include "SmartAssert.h"
#include <cstdio>

/** 
    \brief Create a new graph.
	
    \arg \c N number of nodes.
    \arg \c directional whether the graph should be directional.
*/
Graph::Graph(int N, bool directional)
{
    SMART_ASSERT (N>0)(N);
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

Graph::~Graph()
{
#if GRAPH_DEBUG_LEVEL > 90
    printf("Destroying Graph\n");
#endif
}

/**
   \brief Calculate distances to node j from all nodes.
   
   Uses Dijkstra's algorithm to calculate distances.

   \arg \c dist A pointer to a distance vector for storing the result.
   \arg \c j The target node.
   \arg \c C Constraint vector. If not NULL, then all non-zero entries mean that the algorithm will ignore all corresponding nodes.

   \note To calculate distances \em from node j you only need to
   invert the graph. For undirected graphs, distances from and to any
   node are equal.
*/
bool Graph::CalculateDistances (real* dist, int j, int* C)
{
    for (int n=0; n<N; n++) {
        dist[n] = -1.0;
    }
    dist[j] = 0.0;
    rCalculateDistance_iter(dist, j, C);
    return true;
}

/// Recursive implementation of Dijkstra's algorithm.
bool Graph::rCalculateDistance_iter (real* dist, int j, int* C)
{
    real current_dist = dist[j];
    for (int n=0; n<N; n++) {
        if ((C==NULL)||(C[n]==0)) { 
            if ((n!=j)&&(edge(n,j))) {
                real d = distance(n,j) + current_dist;
                if ((dist[n]<0.0)||(d<dist[n])) {
                    dist[n] = d;
                    //printf ("D[%d]=%f\n", n, d);
                    rCalculateDistance_iter(dist, n, C);
                }
            }
        }
    }
    return true;
}

/// Detecting cylces
bool Graph::hasCycles_iter (std::vector<bool>& mark, int n)
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
bool Graph::hasCycles ()
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

/// Calculate whether j is reachable from all nodes.
bool Graph::isReachable (int j, int* C)
{
    bool flag = true;
    real* D = new real [N];
    CalculateDistances (D, j);
    for (int i=0; i<N; i++) {
        //printf ("D[%d] = %f\n", i, D[i]);
        if (D[i]<0.0) {
            flag = false;
        }
    }

    delete [] D;
    return flag;
}

/// \brief Check whether nodes i and j are k-connected.
///
/// We say that two nodes i and j in an undirected graph are k-connected
/// if there is a path connecting i and j in every subgraph obtained
/// by deleting (k-1) nodes other than i and j together with their adjacent
/// arcs from the graph.
/// The function should always return a number between 0 and N-2.
/// The k-connectivity between a node and itself and that of two nodes
/// connected with a edge is normally undefined. However, here we return
/// the value N-2. 
int Graph::kConnectivity(int i, int j, int* C)
{
    /// \note We assume that two nodes connected with an edge are
    /// (N-2)-connected.  Also, that i and j are k-connected if and only if
    /// either i and j are connected with an arc or there are at least
    /// k node-disjoint paths connecting i and j
    if (i==j) return N-2;
    if (edge(i,j)) {
#if GRAPH_DEBUG_LEVEL >= 99
        printf ("(%d,%d)=edge\n", i,j);
#endif
        return N-2;
    } //i.e. must delete all nodes and at least one of i, j
    int* G;
    G = new int [N];

    if (C) {
        for (int a=0; a<N; a++) {
            G[a] = C[a];
        }
    } else {
        for (int a=0; a<N; a++) {
            G[a] = 0;
        }
    }

    real* dist = new real [N];
    int k = 0;
    if (G[j]) {
        printf ("Err\n");
        delete [] G;
        delete [] dist;
        return 0;
    }
    CalculateDistances (dist, j, G);

    while (MarkShortestPath (dist, G, G, i, j)) {
        k++;
        //printf ("%d\n", k);
    }
    delete [] dist;
    delete [] G;
#if GRAPH_DEBUG_LEVEL >= 99
    printf ("(%d,%d)=%d\n", i,j,k);
#endif
    return k;
}

/// \brief Mark shortest path P, from i to j given constraints and distances.
///
/// Given a distance measure over the set of nodes relative to j
/// mark the shortest path from i to j nodes on P, given constraints C.
/// If there exist many paths of equal length, the first one found is marked.
bool Graph::MarkShortestPath (real* dist, int* C, int* P, int i, int j)
{
    int n = i;
    real min_dist = dist[i];
    int arg_min_dist = i;
    //printf ("[%d]",n);
    for (int m=0; m<N; m++) {
        if (!C[m]) { //don't go to forbidden nodes
            if (edge(n,m)) {
                if (m==j) {
                    return true;
                }
                if (dist[m] < min_dist) {
                    min_dist = dist[m];
                    arg_min_dist = m;
                }
            }
        }
    }
    if (arg_min_dist == n) {
        //printf ("!");
        return false;
    }
    P[arg_min_dist] = 1;
    return MarkShortestPath (dist, C, P, arg_min_dist, j);
}

/// Determine the connectivity of node i
int Graph::kConnectivity(int i, int* C)
{
    int k = N-2;
    for (int j=0; j<N; j++) {
        if ((C==NULL)||(C[j])==0) {
            if (j!=i) {
                int kij = kConnectivity (i, j, C);
                //printf ("kij %d %d=%d\n", i,j, kij);
                if (kij < k) {
                    k = kij;
                }
            }
        }
    }
    if (k>N-2) k = 0;
#if GRAPH_DEBUG_LEVEL >= 99
    printf ("\t<%d>=%d\n", i, k);
#endif
    return k;
}

/// \brief Check whether graph is k-connected.
///
/// Choose an arbitrary node n0 and check k-connectivity between it and
/// every other node. Delete n0 and its adjacent arcs from the graphs,
/// choose another node n1 and check (k-1) connectivity between that
/// node and every other node. Continue in this manner until either
/// node nk-1 is checked to be 1-connected to every remaining node or 
/// (k-i)-connectivity of some node to every remaining nodecannot be verified.
/// In this implementation we use constrain vector C instead of really deleting
/// nodes.
int Graph::kConnectivity()
{
    int* C = new int [N];
    int k = N;
    for (int a=0; a<N; a++) {
        C[a] = 0;
    }
    //Naive implementation
    for (int n=0; n<N; n++) {
        int kn = kConnectivity(n);//,C);
        if (kn < k) {
            k = kn;
        }
    }
#if GRAPH_DEBUG_LEVEL >= 99
    printf ("Naive K:%d\n", k);
#endif
#if 0
    k = N-2;
    for (int n=0; n<N; n++) {
        int kn = kConnectivity(n,C);
        if ((kn==0)||(kn>N-2)) {
            return k-1;
        }
        printf ("(%d)", kn);
        if (kn + n < k) {
            k = kn + n;
        }
        C[n] = 1; //'delete' node
    }
#endif
    delete [] C;
    return k;
}

/// Number of outgoing edges from node i
int Graph::n_out_edges(int i) {
    int cnt = 0;
    for (int j=0; j<N; j++){
        if (edge(i,j)) cnt++;
    }
    return cnt;
}

/// Number of incoming edges to node i
int Graph::n_in_edges(int i) {
    int cnt = 0;
    for (int j=0; j<N; j++){
        if (edge(j,i)) cnt++;
    }
    return cnt;
}

/// Total number of edges for node i
int Graph::n_edges(int i) {
    return n_edges(1);
}


