/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TREE_H
#define TREE_H

#include <list>


struct VoidTreeNode;


struct VoidTreeNodeList
{
    std::list<VoidTreeNode*> next;
    void AddNode(VoidTreeNode* node)
    {
        next.push_back(node);
    }
    ~VoidTreeNodeList()
    {
        next.clear();
    }
};


struct VoidTreeNode
{
    void* data;
    VoidTreeNodeList next;
    VoidTreeNode(void* data_) : data(data_)
    {
    }
};


class VoidTree
{
public:
	VoidTreeNode* root;
	VoidTree(void* node);
	~VoidTree();
};


#endif
