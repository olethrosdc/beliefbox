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

template<class T, class B>
class BeliefExpansionMethod
{
protected:
    T* tree;
    std::vector<T::Node*> node_set;
public:
    virtual ~ExpansionMethod()
    {
    }

    /// Set the reference to a tree
    void setTree(T* tree_)
    {
        tree = tree_;
        node_set = tree->getNodes();
    }

    /// Find the leaf nodes in the node set
    void findLeafNodes(std::vector<int>& leaf_nodes)
    {
        for (uint i=0; i<node_set.size(); ++i) {
            if (node_set[i]->outs.size()==0) {
                leaf_nodes.push_back(i);
            } else {
                n_edge_nodes++;
            }
        }
    }
    virtual int Expand() = 0;
};

class SerialExpansion : public BeliefExpansionMethod
{
    virtual ~SerialExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        assert(leaf_nodes.size() > 0);
        return = leaf_nodes[0];
    }
};


class RandomExpansion : public BeliefExpansionMethod
{
    virtual ~RandomExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        int X =  floor(urandom()*((real) leaf_nodes.size()));
        return leaf_nodes[X];
    }
};


class MeanValueExpansion : public BeliefExpansionMethod
{
    virtual ~MeanValueExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<n_leaf_nodes; i++) {
            BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
            U[i] = node->belief.getGreedyReturn(node->state, gamma);
        }
        return leaf_nodes[ArgMax(U)];
    }
};

class DiscountedMeanValueExpansion : public BeliefExpansionMethod
{
    virtual ~DiscountedMeanValueExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<n_leaf_nodes; i++) {
            BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
            U[i] = ((real) node->depth) * log(gamma)
                + log(node->belief.getGreedyReturn(node->state, gamma));
                            
        }
        return leaf_nodes[ArgMax(U)];
    }
};

class ThompsonSamplingExpansion : public BeliefExpansionMethod
{
    virtual ~ThompsonSamplingExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<n_leaf_nodes; i++) {
            BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
            U[i] = node->belief.sampleReturn(node->state, gamma);
        }
        return leaf_nodes[ArgMax(U)];
    }
};

class DiscountedThompsonSamplingExpansion : public BeliefExpansionMethod
{
    virtual ~DiscountedThompsonSamplingExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        for (int i=0; i<n_leaf_nodes; i++) {
            BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
            U[i] = ((real) node->depth) * log(gamma)
                + log(node->belief.sampleReturn(node->state, gamma));
        }  
        return leaf_nodes[ArgMax(U)];
    }
};


class ThompsonBoundExpansion : public BeliefExpansionMethod
{
    virtual ~ThompsonBoundExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        //std::vector<real> L(leaf_nodes.size());
        for (int i=0; i<n_leaf_nodes; i++) {
            BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
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


class DiscountedThompsonBoundExpansion : public BeliefExpansionMethod
{
    virtual ~ThompsonBoundExpansion()
    {
    }
    virtual int Expand()
    {
        std::vector<int> leaf_nodes;
        findLeafNodes(leaf_nodes);
        std::vector<real> U(leaf_nodes.size());
        std::vector<real> U(leaf_nodes.size());
        //std::vector<real> L(leaf_nodes.size());
        for (int i=0; i<n_leaf_nodes; i++) {
            BeliefTree<SimpleBelief>::Node* node = node_set[leaf_nodes[i]];
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

