/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "CoverTree.h"

/// Constructor needs a point and a level
CoverTree::Node::Node (const CoverTree& tree_, 
                       const Vector& point_,
                       const int level_)
	:  tree(tree_),
       point(point_),
       level(level_), 
       children_level(level_)
{
}

/// Destructor
CoverTree::Node::~Node()
{
	for (uint i=0; i<children.size(); ++i) {
		delete children[i];
	}
}

/// Insert a new point at the given level, as a child of this node
const CoverTree::Node* CoverTree::Node::Insert(const Vector& new_point, const int level)
{

#ifdef DEBUG_COVER_TREE
	printf(" | [%d] ", this->level);
	point.print(stdout);
	printf(" |--(%d)--> ", level);
	new_point.print(stdout);
#endif

	assert(level <= this->level); //hm, does this assert make sense?
	Node* node = new Node(tree, new_point, level);
	children.push_back(node);
	if (level < children_level) {
		children_level = level;
	}
    return node;
}


void CoverTree::Node::Show() const
{
	printf ("%d ", level);
	point.print(stdout);
	if (Size()) {
		printf ("# >>\n");
		for (int i=0; i<Size(); ++i) {
			children[i]->Show();
		}
		printf ("# <<\n");
	} else {
		printf ("# --\n");
	}
}
void CoverTree::Node::Show(FILE* fout) const
{
	printf ("%d ", level);
	point.print(stdout);
	if (Size()) {
		printf ("# >>\n");
		for (int i=0; i<Size(); ++i) {
			children[i]->Show();
		}
		printf ("# <<\n");
	} else {
		printf ("# --\n");
	}
}

const real CoverTree::metric(const CoverSet& Q, const Vector& p) const
{
	real D = INF;
	for (int i=0; i<Q.Size(); ++i) {
		real d_i = metric(Q.nodes[i]->point, p);
		if (d_i < D) {
			D = d_i;
		}
	}
	return D;
}

/** Insert a new point in the tree.
	
	Q_i is the set of points such that the new point may be a nearest
	neighbour to their children.

	If Q_i has depth D, then, for any x, y in Q_i, d(x, y) > 2^D.

	For any x in Q_i, let C(x) be its children. Then,
	for any y in C(x), d(x, y) < 2^D.

	The function is such that Q_i only contains points which whose
	distance to the new point is smaller than 2^level.
*/
const CoverTree::Node* CoverTree::Insert(const Vector& new_point,
                                         const CoverSet& Q_i,
                                         const int level)
{
	Node* closest_node = NULL;
	
	// Check if d(p, Q) > 2^level
	real log_separation = level * log_c;
	real separation = exp(log_separation);
	//Q_i.Show();

	bool separated = true;
	
	// The set of nodes 2^d-close to the new point
	CoverSet Q_next;
	
	// go through all the children and only add them if they are close
    int max_next_level = -INF;
	for (int k=0; k<Q_i.Size(); ++k) {
		int n_children = Q_i.NChildren(k);
		for (int j=-1; j<n_children; ++j) {
			Node* node;
            real dist_i;
			if (j >= 0) {
				node = Q_i.nodes[k]->children[j];
                dist_i = metric(new_point, node->point);
                // ignore children which are too deep.
                if (node->level < level) {
                    max_next_level = std::max(node->level, max_next_level);
                    continue;
                }
			} else {
				node = Q_i.nodes[k];
                dist_i = Q_i.distances[k];

            }
			if (dist_i <= separation) {
				separated = false; 
				Q_next.Insert(node, dist_i);
			}
		}
	}

	// If no points are c^d-close then the point was found previously.
    if (separated) {
        return NULL;
    }

	// Try and see whether the point can be inserted in a subtree
	// Maintain only the points within 2^level distance.
	const Node* found = Insert(new_point, Q_next, max_next_level);

    // The new point x is only possible 
	if (!found) {
		real distance = INF;
		for (int k=0; k<Q_i.Size(); ++k) {
			Node* node = Q_i.nodes[k];
			real dist_k = Q_i.distances[k];//metric(new_point, node->point);
			if (dist_k < distance) {
				distance = dist_k; 
				closest_node = node;
				if (distance <= separation) { // assuming only one node can be here. 
					break;
				}
			}
		}
		
		if (distance <= separation) {
			int new_level = level - 1;
			const Node* inserted = closest_node->Insert(new_point, new_level);
			if (tree_level > new_level) {
				tree_level = new_level;
			}
			return inserted; // Means stop!
		} 
	} 
	return found;
		
}


/// Insert a new point in the tree
const CoverTree::Node* CoverTree::Insert(const Vector& new_point)
{
	if (!root) {
#ifdef DEBUG_COVER_TREE
		printf("Adding root at:");
		new_point.print(stdout);
		printf("\n");
#endif
		root = new Node(*this, new_point, std::numeric_limits<int>::max());
		return root;
	}
	real distance = metric(new_point, root->point);
	int level = 1 + (int) ceil(log(distance) / log_c);
	CoverSet Q;
	Q.Insert(root, distance);
	return Insert(new_point, Q, level);
}


/** Find the nearest node.
   
    If the current node is closest, return that.  
    
    Look through all children which are close enough to this point.
  */
std::pair<const CoverTree::Node*, real> CoverTree::Node::NearestNeighbour(const Vector& query, const real distance) const
{
    std::pair<const CoverTree::Node*, real> retval(this, distance);

	real log_separation = level * tree.log_c;
	real separation = exp(log_separation);

    real& dist = retval.second;
    
    for (int j=0; j<Size(); ++j) {
        real dist_j = children[j]->distanceTo(query);
        if (dist_j - separation <= dist) {
            std::pair<const CoverTree::Node*, real> sub
                = children[j]->NearestNeighbour(query, dist_j);
            //printf ("dist: %f\n", dist_j);
            if (sub.second < dist) {
                retval = sub;
            }
        } else {
            printf("Sep: %f, Dist: %f, Parent: %f, ignoring node [%d -> %d]: ",
                   separation, dist_j, dist, level, children[j]->level);
            children[j]->point.print(stdout);
        }
    }
    //printf("Min dist: %f\n", dist);
	return retval;

}



/// FInd the nearest neighbour in the tree
const CoverTree::Node* CoverTree::NearestNeighbour(const Vector& query_point) const
{
#ifdef DEBUG_COVER_TREE_NN
	printf("Query: "); query_point.print(stdout);
#endif

	if (!root) {
		return NULL;
	}

	std::pair<const CoverTree::Node*, real> val
        = root->NearestNeighbour(query_point, root->distanceTo(query_point));
    return val.first;
}

/** Check that the tree implements the constraints properly */
bool CoverTree::Check() const
{
	if (!root) {
		return true;
	}
	CoverSet Q;
	Q.Insert(root, INF);
	return Check(Q, std::numeric_limits<int>::max());
}

/** Internal check method.

	If parents are at a particular level, then it is necessary
	for all of them to have a separation s > 2^d.
*/
bool CoverTree::Check(const CoverSet& parents, const int level) const
{
	real separation = Separation(parents);
	if (log(separation - 2) <= static_cast<real>(level)) {
		return false;
	}
	return true;
}

real CoverTree::Separation(const CoverSet& Q) const
{
	real separation = INF;
	for (int i=0; i<Q.Size(); ++i) {
		for (int j=i+1; j<Q.Size(); ++j) {
			real s_ij = metric(Q.nodes[i]->point, Q.nodes[j]->point);
			if (s_ij < separation) {
				separation = s_ij;
			}
		}
	}
	return separation;
}

/** Show the tree
 */
void CoverTree::Show() const
{
	if (root) {
		FILE* fout = fopen("tree.dot", "w");
		if (!fout) {
			fprintf(stderr, "Could not write to dot file\n");
			exit(-1);
		}
		fprintf (fout, "digraph Covertree {\n");
		fprintf (fout, "ranksep=2; rankdir=TB; \n");
		root->Show();
		fprintf (fout, "}\n");
		fclose(fout);
	} else {
		printf ("# Tree is empty\n");
	}
}

/** Default constructor */
CoverTree::CoverTree(real c)
{
    assert (c > 0);
    log_c = c;
	root = NULL;
	tree_level = std::numeric_limits<int>::max();
}

/** Destructor */
CoverTree::~CoverTree()
{
	delete root;
}

