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


//CoverTree::Statistics::Statistics(const CoverTree& tree_, const real& w, int i_dim, int o_dim, real N0, real a, Matrix S0) 
//		: tree(tree_) 
//{
//	StatePrediction		= new BayesianMultivariateRegression(i_dim, o_dim, (N0 * Matrix::Unity(o_dim, o_dim)), N0, a);
//	RewardPrediction	= new BayesianMultivariateRegression(i_dim, 1, (N0 * Matrix::Unity(1, 1)), N0, a);
//	weight = w;
//}
//
///// Deconstructor
//CoverTree::Statistics::~Statistics()
//{
//	delete StatePrediction;
//	delete RewardPrediction;
//}

/// Constructor needs a point and a level
CoverTree::Node::Node (const CoverTree& tree_, 
					   const Vector& point_,
                       const int level_, 
					   CoverTree::Node* const father_,
					   void* object_)
	:	tree(tree_),
		point(point_),
		level(level_), 
		children_level(level_),
		father(father_),
		object(object_)
{	
//	point.print(stdout);
	index = tree.GetNumNodes();
	descendants		= 0;
	samples			= 0;
	active_flag		= false;
	active_index	= -1;
	basis_flag		= false;
	basis_index		= -1;
	
	child			= NULL;

	weight = 1;
	depth  = 1;
	if(father!=NULL) {
		weight	= father->GetWeight() / 2.0;
		depth	= father->Depth() + 1;
	}
	
	int i_dim;
	if(tree.RBFs) {
		i_dim = tree.RBFs->size() + 1;
	}
	else {
		i_dim = point.Size() + 1;
	}
	int o_dim = point.Size();
	real N0		= tree.N0;
	real a		= tree.a;

	StatePrediction		= new BayesianMultivariateRegression(i_dim, o_dim, (N0*Matrix::Unity(o_dim, o_dim)), N0, a);
	if(tree.RewardPred) {
		RewardPrediction = new BayesianMultivariateRegression(i_dim, 1, (N0*Matrix::Unity(o_dim, o_dim)), N0, a);
	}
}

/// Destructor
CoverTree::Node::~Node()
{
	for (uint i=0; i<children.size(); ++i) {
		delete children[i];
	}

	delete StatePrediction;
	if(tree.RewardPred) 
		delete RewardPrediction;
}

/// Insert a new point at the given level, as a child of this node
void CoverTree::Node::Insert(const Vector& new_point, const int level, void* obj)
{

#ifdef DEBUG_COVER_TREE
	printf(" | [%d] ", this->level);
	point.print(stdout);
	printf(" |--(%d)--> ", level);
	new_point.print(stdout);
#endif

	assert(level <= this->level); //hm, does this assert make sense?
	Node* node	 = new Node(tree, new_point, level, this, obj);
//	node->father = this;
	children.push_back(node);
	if (level < children_level) {
		children_level = level;
	}
}

/// Insert a new point at the given level, as a child of this node
void CoverTree::Node::Insert(const Vector& new_point, const Vector& phi, const Vector& next_state, const real& reward, const bool& absorb, const int level, void* obj)
{
#ifdef DEBUG_COVER_TREE
	printf(" | [%d] ", this->level);
	point.print(stdout);
	printf(" |--(%d)--> ", level);
	new_point.print(stdout);
#endif
	
	assert(level <= this->level); //hm, does this assert make sense?
	Node* node = new Node(tree, new_point, level, this, obj);
//	node->father = this;
	children.push_back(node);
	node->UpdateStatistics(phi, next_state, reward, absorb, NULL);
	if (level < children_level) {
		children_level = level;
	}
}

void CoverTree::Node::Show() const
{
	printf ("level: %d ", level);
	printf("descendants: %d ", descendants);
	printf("depth: %d ", depth);
	printf("index: %d ", index);
	printf("#children: %d", Size());
	printf("Weights: %f ", weight );
	printf("Obs: %d ", NumObs());
	printf("Active: %s ",(active_flag)?"true":"false");
	printf("Active_index: %d ", active_index);

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
void CoverTree::Node::ShowSampling() const
{
	if(active_flag == true) {
//		printf ("level: %d ", level);
		point.print(stdout);
		printf ("level: %d ", level);
		printf("descendants: %d ", descendants);
		printf("depth: %d ", depth);
		printf("samples: %d ", samples);
		printf("index: %d ", index);
		if(active_flag == true) {
			if (Size()) {
//				printf ("# >>\n");
				for (int i=0; i<Size(); ++i) {
					children[i]->ShowSampling();
				}
//				printf ("# <<\n");
			} else {
//				printf ("# --\n");
			}
		}
	}
}
void CoverTree::Node::ShowBasis() const
{
	if(GetActiveBasis() == true) {
		point.print(stdout);
		if (Size()) {
			//			printf ("# >>\n");
			for (int i=0; i<Size(); ++i) {
				children[i]->ShowBasis();
			}
			//			printf ("# <<\n");
		} //else {
//			//			printf ("# --\n");
//		}
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
void CoverTree::Node::ShowSampling(FILE* fout) const
{
	printf ("level: %d ", level);
	point.print(stdout);
	if(active_flag) {
		if (Size()) {
			printf ("# >>\n");
			for (int i=0; i<Size(); ++i) {
				children[i]->ShowSampling();
			}
			printf ("# <<\n");
		} else {
			printf ("# --\n");
		}
	}
}
void CoverTree::Node::ShowBasis(FILE* fout) const
{
	printf ("level: %d ", level);
	point.print(stdout);
	if(basis_flag) {
		if (Size()) {
//			printf ("# >>\n");
			for (int i=0; i<Size(); ++i) {
				children[i]->ShowBasis();
			}
//			printf ("# <<\n");
		} else {
//			printf ("# --\n");
		}
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

void CoverTree::UpdateStatistics(const Vector& input, const Vector& output) {
	real total_probability = 0;
//	printf("ROOT\n");
	Node* examine = root; 
	while(examine != NULL) {
		total_probability = examine->Update(input, output, total_probability);
		examine = examine->child;
	}
}

const void CoverTree::SamplingNode(Node* n) 
{
	///In this point, we sample the regressor's nodes.
	if(n == root) {
		num_sampling_nodes++;
		n->active_flag	= true;
		n->active_index	= GetNumSamplingNodes();
	}
	else {
		if(urandom() <= n->weight && n->samples > thres) {
//		printf("Samples = %d, Depth = %d\n",n->samples, n->depth);
//		Vector nn = n->point;
//		nn.print(stdout);
			num_sampling_nodes++;
			n->active_flag	= true;
			n->active_index	= GetNumSamplingNodes();
		}
		else {
			n->active_flag	= false;
			n->active_index	= -1;
		}
		///In this point, we sample the nodes which play the basis function role.
		if(n->samples > GetBasisThreshold()) {
			num_basis_nodes++;
			n->basis_flag = true;
			n->basis_index = num_basis_nodes;
			Basis.push_back(n);
		}
		else {
			n->basis_flag  = false;
			n->basis_index = -1;
		}
	}
	if(n->Size()) {
		for(int i=0; i<n->Size(); ++i) {
			SamplingNode(n->children[i]);
		}
	}
}

const void CoverTree::SamplingNode(Node* n, const Vector& R)
{
	if(R[n->index] <= n->weight && n->samples > thres) {
		num_sampling_nodes++;
		n->active_flag	= true;
		n->active_index	= GetNumSamplingNodes();
	}
	else {
		n->active_flag		= false;
		n->active_index		= -1;
	}
	///In this point, we sample the nodes which play the basis function role.
	if(n->samples > GetBasisThreshold()) {
		num_basis_nodes++;
		n->basis_flag = true;
		n->basis_index = num_basis_nodes;
		Basis.push_back(n);
	}
	else {
		n->basis_flag  = false;
		n->basis_index = -1;
	}
	
	if(n->Size()) {
		for(int i=0; i<n->Size(); ++i) {
			SamplingNode(n->children[i],R);
		}
	}
}

const void CoverTree::SamplingTree()
{
	SamplingNode(root);
}

const void CoverTree::SamplingTree(const Vector& R)
{
	assert(R.Size() == num_nodes);
	SamplingNode(root,R);
}

Vector CoverTree::BasisCreation(const Vector& state) const
{
	Vector phi = state;
	if(RBFs != NULL) {
		RBFs->Evaluate(state);
		phi = RBFs->F();
	}
	int dim = phi.Size() + 1;
	phi.Resize(dim);
	phi[dim-1] = 1.0;
	return phi;
}

const std::vector< std::pair<int,real> > CoverTree::ExternalBasisCreation(const Vector& state) const 
{
	std::vector< std::pair<int, real> > phi;
	std::pair<int, real> path_data;
	const Node *path_node = NearestNeighbour(state);
//	real total_weights = 1;
	
	while(path_node != NULL) {
		if(path_node->GetActiveBasis() == true) {
//		if(path_node->GetActive() == true) {
//			printf("adasdfsfasddsaffffffff\n");
//     		path_data.first	 = path_node->GetSamplingIndex();
//			path_data.second = path_node->GetWeight() * total_weights;
//		    path_data.first  = path_node->GetIndex();
//			path_data.second = 1.0;
//			printf("Index = %d  === Weight = %f\n",path_node->GetSamplingIndex(), path_data.second);
			path_data.first	 = path_node->GetBasisIndex();
			real beta = pow(2.0, (real)path_node->level);
			Vector d = pow((state - path_node->point)/beta, 2.0);
			real r = d.Sum();
			//real d = EuclideanNorm(&x, &center);
			path_data.second = exp(-0.5*r);
//			path_data.second = 1.0;
			phi.push_back(path_data);
//			total_weights *= (1 - path_node->GetWeight());
//			printf("adasdfsfasddsaffffffff\n");
		}
		path_node = path_node->father;
	}
//	for(uint i = 0; i<Basis.size(); ++i) {
//		real beta = pow(2.0, (real)Basis[i]->level) / 1.0;
//		Vector d = pow((state - Basis[i]->point)/beta, 2.0);
//		real r = d.Sum();
//		path_data.first = Basis[i]->GetBasisIndex();
//		path_data.second = exp(-0.5*r);
//		phi.push_back(path_data);
//	}
	return phi;
}

/** Insert a new point in the tree.
	
	Q_i is the set of points such that the new point may be a nearest
	neighbour to their children.

	If Q_i has depth D, then, for any x, y in Q_i, d(x, y) > 2^D.

	For any x in Q_i, let C(x) be its children. Then,
	for any y in C(x), d(x, y) < 2^D.

	The function is such that Q_i only contains points which whose
	distance to the new point is smaller than 2^level.

    Officially:

    Insret (p, Q, i)
    Q = {C(q) : q \in Q_i}
    if (d(p, Q) > 2^i) {
       return parent found (true)
    } else  {
       Q_{i-1} = {q in Q: d(p, q) \leq 2^i}
       found = Insert(p, Q_{i-1}, i - 1)
       if (found and d(p, Q_i) \leq 2^i) {
          pick a single (any) q in Q_i s.t. d(p,q) \leq 2^i
          insert p in C(q)
          return finished (false)
       } else {
          return found;
       }
   }
*/
CoverTree::Node* CoverTree::Insert(const Vector& new_point,
                                         const CoverSet& Q_i,
                                         const int level,
										 void* obj)
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
            real dist_i = 0;
			if (j >= 0) {
				node = Q_i.nodes[k]->children[j];
                // ignore children which are too deep.
                if (node->level < level) {
                    max_next_level = std::max(node->level, max_next_level);
                    continue;
                }
                dist_i = metric(new_point, node->point);
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
	Node* found = Insert(new_point, Q_next, max_next_level, obj);

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
 		
		if (distance <= separation && GetEntranceThreshold(closest_node->Depth()) <= closest_node->NumObs()) {
			int new_level = level - 1;
			num_nodes++;
			closest_node->Insert(new_point, new_level);
			closest_node->descendants = closest_node->descendants + 1;
			if (tree_level > new_level) {
				tree_level = new_level;
			}
			return closest_node; //Means stop!!
		} 
		else if(distance <= separation && GetEntranceThreshold(closest_node->Depth()) >= closest_node->NumObs()) {
			return closest_node; // Means stop!!!!!!!!!!
		}
	} 
	else if(found->father != NULL){
		if(found->father != found){
			if(found->descendant_flag) {
				found->descendant_flag = false;
				found = found->father;
				found->descendants++;
				found->descendant_flag = true;
			}
			else {
				found = found->father;
			}
		}
	}
	return found;
		
}

CoverTree::Node* CoverTree::Insert(const Vector& new_point, const Vector& phi, const Vector& next_state, const real& reward, const bool& absorb, const CoverSet& Q_i, const int level, void* obj)
{
	Node* closest_node = NULL;
	
	// Check if d(p, Q) > 2^level
	real log_separation = level * log_c;
	real separation = exp(log_separation);
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
                // ignore children which are too deep.
                if (node->level < level) {
                    max_next_level = std::max(node->level, max_next_level);
                    continue;
                }
                dist_i = metric(new_point, node->point);
			} else {
				node = Q_i.nodes[k];
                dist_i = Q_i.distances[k];
            }
			if (dist_i <= 1e-3) {// CUTOFF 
				return node;
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
	Node* found = Insert(new_point, phi, next_state, reward, absorb, Q_next, max_next_level, obj);
	
    // The new point x is only possible 
	if (!found) {
		real distance = INF;
		for (int k=0; k<Q_i.Size(); ++k) {
			Node* node = Q_i.nodes[k];
			real dist_k = Q_i.distances[k];    //metric(new_point, node->point);
			if (dist_k < distance) {
				distance = dist_k; 
				closest_node = node;
				if (distance <= separation) { // assuming only one node can be here. 
					break;
				}
			}
		}
		
		if (distance <= separation && GetEntranceThreshold(closest_node->Depth()) <= closest_node->NumObs()) {
			int new_level = level - 1;
			num_nodes++;
			closest_node->Insert(new_point, phi, next_state, reward, absorb, new_level);
			closest_node->descendants++;
			closest_node->descendant_flag = true;
			if (tree_level > new_level) {
				tree_level = new_level;
			}
			return closest_node; // Means stop!
		} 
		else if(distance <= separation && GetEntranceThreshold(closest_node->Depth()) >= closest_node->NumObs()) {
			closest_node->UpdateStatistics(phi, next_state, reward, absorb, NULL);
			return closest_node; // Means stop!
		}
	}
	else if(found->father != NULL) {
		if(found->father != found){
			Node* child = found;
			if(found->descendant_flag) {
				found->descendant_flag = false;
				found = found->father;
				found->descendants++;
				found->descendant_flag = true;
				found->UpdateStatistics(phi, next_state, reward, absorb, child);
			}
			else {
				found = found->father;
				found->UpdateStatistics(phi, next_state, reward, absorb, child);
			}
		}
	}
	return found;		
}

/// Insert a new point in the tree
CoverTree::Node* CoverTree::Insert(const Vector& new_point, void* obj)
{
	total_samples++;
	if (!root) {
#ifdef DEBUG_COVER_TREE
		printf("Adding root at:");
		new_point.print(stdout);
		printf("\n");
#endif
		num_nodes = 1;
		root = new Node(*this, new_point, std::numeric_limits<int>::max(), NULL, obj);
		return root;
	}
	real distance = metric(new_point, root->point);
	int level = 1 + (int) ceil(log(distance) / log_c);
	CoverSet Q;
	Q.Insert(root, distance);
	return Insert(new_point, Q, level, obj);
}

/// Insert a new point in the tree along with its statistics
CoverTree::Node* CoverTree::Insert(const Vector& new_point, const Vector& next_state, const real& reward, const bool& absorb, void* obj)
{
	Vector phi = BasisCreation(new_point);
	total_samples++;
	if (!root) {
#ifdef DEBUG_COVER_TREE
		printf("Adding root at:");
		new_point.print(stdout);
		printf("\n");
#endif
		num_nodes = 1;
		root = new Node(*this, new_point, std::numeric_limits<int>::max(), NULL, obj);

		root->UpdateStatistics(phi, next_state, reward, absorb, NULL);
		return root;
	}
	
	real distance = metric(new_point, root->point);
	int level = 1 + (int) ceil(log(distance) / log_c);
	
	CoverSet Q;
	Q.Insert(root, distance);
	
	///A new point is inserted in our tree
	CoverTree::Node* inserted =  Insert(new_point, phi, next_state, reward, absorb, Q, level, obj);

	///Update path statistics.
	UpdateStatistics(phi, next_state);

	return inserted;
}

/** Find the nearest node.
   
    If the current node is closest, return that.  
    
    Look through all children which are close enough to this point.
  */
std::pair<const CoverTree::Node*, real> CoverTree::Node::NearestNeighbour(const Vector& query, const real distance) const
{
    std::pair<const CoverTree::Node*, real> retval(this, distance);

	//real log_separation = level * tree.log_c;
	//real separation = exp(log_separation);

    real& dist = retval.second;
    
    for (int j=0; j<Size(); ++j) {
        real dist_j = children[j]->distanceTo(query);
        //        if (dist_j - separation <= dist) {
        if (dist_j <= dist + exp(children[j]->level * tree.log_c)) {
            std::pair<const CoverTree::Node*, real> sub
                = children[j]->NearestNeighbour(query, dist_j);
//            printf ("dist: %f\n", dist_j);
            if (sub.second < dist) {
                retval = sub;
            }
            //} else {
            //            printf("Sep: %f, Dist: %f, Parent: %f, ignoring node [%d -> %d]: ",
            //separation, dist_j, dist, level, children[j]->level);
            //children[j]->point.print(stdout);
        }
    }
//    printf("Min dist: %f\n", dist);
	return retval;

}

/// Find the nearest neighbour in the tree
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

const CoverTree::Node* CoverTree::SelectedNearestNeighbour(const Vector& query_point) const
{
	const Node* found = NearestNeighbour(query_point);
	while(found!=NULL) {
		if(found->active_flag) {
			return(found);
		}
		found = found->father;
	}
	return found;
}
const Vector CoverTree::GenerateState(const Vector& query_point) const
{
	const Node* found	= SelectedNearestNeighbour(query_point);
	Vector phi			= BasisCreation(query_point);
	return found->GenerateS(phi);
}
const real CoverTree::GenerateReward(const Vector& query_point) const
{
	const Node* found = SelectedNearestNeighbour(query_point);
	Vector phi			= BasisCreation(query_point);
	return found->GenerateR(phi);
}
const int CoverTree::GetNumNodes() const {
	return num_nodes;
}

const void CoverTree::SetNumNodes(int num) {
	num_nodes = num;
}

const int CoverTree::GetNumSamplingNodes() const {
	return num_sampling_nodes;
}

const void CoverTree::SetNumSamplingNodes(int num) {
	num_sampling_nodes = num;
}

const int CoverTree::GetNumBasisNodes() const {
	return num_basis_nodes;
}

const void CoverTree::SetNumBasisNodes(int num) {
	num_basis_nodes = num;
}

const void CoverTree::Reset() {
	delete root;
	root		= NULL;
	num_nodes			= 0;
	num_sampling_nodes	= 0;
	num_basis_nodes		= 0;
	total_samples		= 0;
	num_nodes			= 0;
	Basis.clear();
	tree_level	= std::numeric_limits<int>::max();
}

const real CoverTree::GetEntranceThreshold(int depth) const {
	// Make sure we have enough observations to justify adding a new node. 
	// This means at least as many as total outcomes.
    // real threshold = (real) n_outcomes; 
	
    // Go deeper when there has been at least one observations node. 
//    real threshold = sqrt((real) depth);//2;
//	real threshold = log((real) depth) + 1;//2;
//	real threshold = pow((real) depth, 2.0);
	real threshold = 0.0;
	
	return threshold;
}

const real CoverTree::GetBasisThreshold() const {
	real n = 1.0 / 10.0;
	return pow((real) total_samples, n); 
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

/** Show the tree
 */
void CoverTree::ShowSampling() const
{
	if (root) {
		FILE* fout = fopen("tree.dot", "w");
		if (!fout) {
			fprintf(stderr, "Could not write to dot file\n");
			exit(-1);
		}
		fprintf (fout, "digraph Covertree {\n");
		fprintf (fout, "ranksep=2; rankdir=TB; \n");
		root->ShowSampling();
		fprintf (fout, "}\n");
		fclose(fout);
	} else {
		printf ("# Tree is empty\n");
	}
}

void CoverTree::ShowBasis() const 
{
	if (root) {
		FILE* fout = fopen("tree.dot", "w");
		if (!fout) {
			fprintf(stderr, "Could not write to dot file\n");
			exit(-1);
		}
		fprintf (fout, "digraph Covertree {\n");
		fprintf (fout, "ranksep=2; rankdir=TB; \n");
		root->ShowBasis();
		fprintf (fout, "}\n");
		fclose(fout);
	} else {
		printf ("# Tree is empty\n");
	}
}
/** Default constructor */
CoverTree::CoverTree(real c, real a_, real N0_, RBFBasisSet* RBFs_, bool f)
{
	assert (c > 1);
    log_c		= log(c);
	RewardPred	= f;
	a			= a_;
	N0			= N0_;
	RBFs		= RBFs_;
	root		= NULL;
	num_nodes			= 0;
	num_sampling_nodes	= 0;
	num_basis_nodes		= 0;
	total_samples		= 0;
	thres		= 200;
	Basis.clear();
	tree_level	= std::numeric_limits<int>::max();
}

/** Destructor */
CoverTree::~CoverTree()
{
	delete root;
}

