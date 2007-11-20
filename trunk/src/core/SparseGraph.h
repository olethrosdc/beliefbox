/* -*- Mode: C++; -*- */
// copyright (c) 2006-2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef SPARSE_GRAPH_H
#define SPARSE_GRAPH_H
#include "Graph.h"
#include <vector>
#include <list>


class SparseGraph : public Graph
{
protected:
    std::vector<HalfEdgeList> parents;
    std::vector<HalfEdgeList> children;
    HalfEdge& getParentHalfEdge(int src, int dst);
    HalfEdge& getChildHalfEdge(int src, int dst);
    virtual bool rCalculateDistance_iter(real* dist, int j, int* C=NULL);
    virtual bool hasCycles_iter (std::vector<bool>& mark, int n);
public:
	SparseGraph(int N, bool directional);
	virtual ~SparseGraph() {}
    bool AddEdge(Edge e, bool clear = false);
	virtual bool edge (int src, int dst); ///< Returns true if there is an edge from \c src to \c dst.
	virtual real distance (int src, int dst); ///< Distance between two nodes. Note that in some applications this would have the meaning of an edge weight rather than a distance.
    virtual HalfEdgeListIterator getFirstParent(int node);
    virtual HalfEdgeListIterator getFirstChild(int node);
    virtual int n_parents(int node);
    virtual int n_children(int node);
};


//class GeneralGraph : public SparseGraph
//{
//   
//};

#endif
