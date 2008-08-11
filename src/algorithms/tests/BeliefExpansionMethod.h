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

#ifndef EXPANSION_METHOD_H
#define BELIEF_EXPANSION_METHOD_H

#include "bandit_uct.h"
#include <vector>
#include <cmath>

template<class T, class Node, class B>
class BeliefExpansionMethod 
{
protected:
    T* tree;
    std::vector<Node> node_set;
    //std::vector<T**> node_set;
public:
    virtual ~BeliefExpansionMethod()
    {
    }

    /// Set the reference to a tree
    void setTree(T* tree_)
    {
        tree = tree_;
        node_set = tree->getNodes();
    }

    virtual void findLeafNodes(std::vector<int>& leaf_nodes)
    {
        for (uint i=0; i<node_set.size(); ++i) {
            if (node_set[i]->outs.size()==0) {
                leaf_nodes.push_back(i);
            }
        }
    }

    virtual int Expand() = 0;
};

template <class T, class Node, class B>
class SerialExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~SerialExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->findLeafNodes(leaf_nodes);
        assert(leaf_nodes.size() > 0);
        return leaf_nodes[0];
    }
};


template <class T, class Node, class B>
class RandomExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~RandomExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        int X =  floor(urandom()*((real) leaf_nodes.size()));
        return leaf_nodes[X];
    }
};


template <class T, class Node, class B>
class MeanValueExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~MeanValueExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<leaf_nodes.size(); i++) {
            Node node = this->node_set[leaf_nodes[i]];
            U[i] = node->belief.getGreedyReturn(node->state, gamma);
        }
        return leaf_nodes[ArgMax(U)];
    }
};

template <class T, class Node, class B>
class DiscountedMeanValueExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~DiscountedMeanValueExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<leaf_nodes.size(); i++) {
            Node node = this->node_set[leaf_nodes[i]];
            U[i] = ((real) node->depth) * log(gamma)
                + log(node->belief.getGreedyReturn(node->state, gamma));
                            
        }
        return leaf_nodes[ArgMax(U)];
    }
};

template <class T, class Node, class B>
class ThompsonSamplingExpasnion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~ThompsonSamplingExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<leaf_nodes.size(); i++) {
            Node* node = node_set[leaf_nodes[i]];
            U[i] = node->belief.sampleReturn(node->state, gamma);
        }
        return leaf_nodes[ArgMax(U)];
    }
};

template <class T, class Node, class B>
class DiscountedThompsonSamplingExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~DiscountedThompsonSamplingExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<leaf_nodes.size(); i++) {
            Node* node = node_set[leaf_nodes[i]];
            U[i] = ((real) node->depth) * log(gamma)
                + log(node->belief.sampleReturn(node->state, gamma));
        }  
        return leaf_nodes[ArgMax(U)];
    }
};

template <class T, class Node, class B>
class ThompsonBoundExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~ThompsonBoundExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        //std::vector<real> L(leaf_nodes.size());
        for (int i=0; i<leaf_nodes.size(); i++) {
            Node* node = node_set[leaf_nodes[i]];
            real Ui = node->belief.sampleReturn(node->state, gamma);
            real Li = node->belief.getGreedyReturn(node->state, gamma);
            if (Ui < Li) {
                Ui = Li;
            }
            U[i] = Ui;
                //L[i] = Li;
        } 
        return leaf_nodes[ArgMax(U)];
    }
};

template <class T, class Node, class B>
class DiscountedThompsonBoundExpansion : public BeliefExpansionMethod<T, Node, B>
{
    virtual ~ThompsonBoundExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        this->LeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        std::vector<real> U(leaf_nodes.size());
        //std::vector<real> L(leaf_nodes.size());
        for (int i=0; i<leaf_nodes.size(); i++) {
            Node* node = node_set[leaf_nodes[i]];
                real Ui = node->belief.sampleReturn(node->state, gamma);
                real Li = node->belief.getGreedyReturn(node->state, gamma);
                if (Ui < Li) {
                    Ui = Li;
                }
                U[i] =  ((real) node->depth) * log(gamma) + log(Ui);
                //L[i] = Li;
        }
        return leaf_nodes[ArgMax(U)];
    }
    
};



#endif

