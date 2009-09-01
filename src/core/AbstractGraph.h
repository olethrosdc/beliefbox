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

#ifndef ABSTRACT_GRAPH_H
#define ABSTRACT_GRAPH_H

#include <cstdlib>
#include <list>
#include <vector>
#include <set>
#include "real.h"

/*@{*/
/// A general graph class.  This class provides an interface for
/// classes implementing graph semantics. Most of the non-pure virtual
/// functions do not need to be over-ridden, unless a sub-class has a
/// special fast implementation for one of them.
class AbstractGraph
{
protected:
    int N; ///< number of nodes
    bool directional; ///< whether it is directional
    virtual bool hasCycles_iter (std::vector<bool>& mark, int n);
public:
    AbstractGraph(int N, bool directional);
    virtual ~AbstractGraph();
    
    /// The number of nodes
    int n_nodes()
    {
        return N;
    }

    /// Whether the graph is directional
    bool is_directional()
    {
        return directional;
    }
    
    virtual int n_out_edges (int i);
    virtual int n_in_edges (int i);
    virtual int n_edges (int i);
    virtual int kConnectivity(); 
    virtual int kConnectivity(int i, int* C = NULL); 
    virtual int kConnectivity(int i, int j, int* C = NULL); 
    virtual bool hasCycles ();
    virtual bool edge (int src, int dst) = 0;
};


struct StructuralNode;

struct StructuralEdge
{
    StructuralNode* src;
    StructuralNode* dst;
};

struct StructuralNode
{
    std::list<StructuralEdge*> out_edges;
    std::list<StructuralEdge*> in_edges;
};



template <typename V, typename Q>
class TemplatedGraph : public AbstractGraph
{
public:
    /// The templated node class
    struct Node
    {
        StructuralNode node;
        V v;
    };

    /// Te templated edge class
    struct Edge
    {
        StructuralEdge edge;
        Q q;
    };
public:
    /// Returns true if there is an edge from \c src to \c dst.
    virtual bool edge (int src, int dst) = 0; 

    /// Value of a node
    virtual V Value (int node) = 0; 

    /// Value of an edge
    virtual Q Value (int src, int dst) = 0;
};

/*@}*/

#endif
