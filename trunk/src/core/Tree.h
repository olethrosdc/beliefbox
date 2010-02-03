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

class VoidTree::Node;

class VoidTree
{
public:
	struct NodeList
	{
		std::list<Node*> next;
		AddNode(Node* node)
		{
			next.push_back(node);
		}
	};
	struct Node
	{
		
	};
	void* root;
	VoidTree();
	~VoidTree
};

template <typename T>
class Tree<T>  
{
	
};

#endif
