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

#ifndef CONNECTIVITY_MATRIX_H
#define CONNECTIVITY_MATRIX_H

#include "Graph.h"
#include "SparseGraph.h"

/// A graph implemented as a connectivity matrix.
class ConnectivityMatrix : public Graph
{
protected:
    int** C; ///< connection matrix
    real** D; ///< distance matrix
public:
    /** \brief Create a matrix with N nodes, directional, with an optional list of prior connections and optional inversion of directions.
     */
    ConnectivityMatrix(int N, bool directional, int** c = NULL, bool invert=false, real** d = NULL);
    ConnectivityMatrix(SparseGraph& G);
    ~ConnectivityMatrix();
    virtual bool edge (int src, int dst);
    virtual real distance (int src, int dst);
    virtual bool CheckSymmetry (); ///< Check whether graph is symmetric
};

#endif
