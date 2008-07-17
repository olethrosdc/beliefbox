// -*- Mode: C++; -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "PolicyEvaluation.h"
#include "BetaDistribution.h"
#include "Random.h"

#include <list>


template <typename B>
class BeliefTree
{
public:
    class Edge;
    class Node
    {
    public:
        B belief;
        int state;
        std::list<Edge&> edges;
    };
    
    class Edge
    {
    public:
        Node& src;
        Node& dst;
    };

    std::vector<Node> N;
    std::vector<Edge> E;
    /// Return 
    int Expand(int i, int a, real r, int s)
    {
    }
    
};



/// A toy UCT stopping problem




//void EvaluateAlgorithm(BeliefExpansionAlgorithm& algorithm, real mean_r);


int main (int argc, char** argv)
{
    real alpha = 0;
    real beta = 0;
    //ConjugatePrior* prior = new BetaDistribution(alpha, beta);

    real mean_r = urandom(-1, 1);
    
    //EvaluateAlgorithm(algorithm, mean_r);
    
    return 0;
}


#endif
