// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BACKWARDS_INDUCTION_H
#define BACKWARDS_INDUCTION_H


#include "SparseGraph.h"
#include "Vector.h"
#include <stdexcept>
#include <iostream>

class InductionOperator
{
public:
    virtual real Induct(real x, real y)= 0;
    virtual ~InductionOperator() {}
};

class MaxInduction : public InductionOperator
{
public:
    virtual real Induct(real x, real y)
    {
	return std::max(x, y);
    }
};

class AddInduction : public InductionOperator
{
public:
    virtual real Induct(real x, real y)
    {
	return x + y;
    }
};


/** The backwards induction algorithm.  It needs a sort of Graph.
    This precludes it from working in anything other than a
    discrete space.
*/
//template <typename StateType, typename ActionType>
template <class VectorType, class GraphType, class InductionType>
class BackwardsInduction
{
protected:
    VectorType& V; ///< A vector where to hold the values
    GraphType& G; ///< The graph on which to do it
    InductionType& J; ///< The type of induction to perform
public:
    BackwardsInduction(VectorType& V_, GraphType& G_, InductionType& J_) :
	V(V_), G(G_), J(J_)
    {}

    void calculate(NodeSet T) {
	if (G.hasCycles()) {
		throw std::domain_error("Graph has loops, cannot use run\n");
	    }
    
	if (T.size() == 0) {
	    std::cerr << "Warning: No terminal nodes specified!" << std::endl;
	}

	while(T.size()) {
	    NodeSet P;
	    for (NodeSetIterator i=T.begin(); i!=T.end(); ++i) {
		int n_parents = G.n_parents(*i);
		int parent = 0;
		for (HalfEdgeListIterator e= G.getFirstParent(*i);
		     parent < n_parents;
		     ++e, parent++) {
		    int j = e->node;
		    //P.push_back(j);
		    P.insert(j);
		    V[j] = J.Induct(e->w, V[i]);
		}
	    }
	    T = P;
	}
    }
    void calculate_loopy(NodeList T, int max_iter, real end_accuracy);
};



#endif

