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

#ifndef COVER_TREE_H
#define COVER_TREE_H

#include "real.h"
#include "Vector.h"
#include "Matrix.h"
#include "Random.h"
#include "BasisSet.h"
#include "BayesianMultivariateRegression.h"

#undef DEBUG_COVER_TREE
#undef DEBUG_COVER_TREE_NN

/** A cover tree.
 
 Implements the algorithm "Cover trees for Nearest Neighbor",
 Beygelzimer, Kakade, Langford. Based on the technical report "Fast
 Nearest Neighbors" of Thomas Kollar.
 
 The tree is constructed from a set \f$S\f$ of points, a metric
 \f$\rho\f$ and a constant \f$c > 1\f$. We write \f$d_{i,j} =
 \rho(x_i, x_j)\f$ for the distance between the i-th and j-th
 points. Let \f$S_n\f$ denote the set of nodes at level \f$n\f$ of
 the tree and \f$C(i)\f$ the set of children of the i-th node. The
 tree has the following properties.
 
 1. If \f$i \in S_n\f$ then \f$i \in S_{n-1}\f$.
 
 2. \f$S_\infty = \fS$.
 
 3. For any \f$i, j \in S_k\f$, \f$d_{i,j} > c^d\f$.
 
 4. If \f$i \in S_n\f$ and \f$j \in C(i)\f$, 
 
 */
class CoverTree
	{
	public:
		/// The raw metric used between points.
		static const real metric(const Vector& x, const Vector& y)
		{
			return L1Norm(x, y);
		}
		/// This simply is a node
		struct Node
		{
			const CoverTree& tree;
			
			Vector point;	       	///< The location of the point.
			std::vector<Node*> children;///< Pointer to children
			int level;;	                ///< Level in the tree
			int index;	       		///< Node index
			int depth;		       	///< Depth in the tree
			int descendants;        ///< Number of descendants
			int samples;	       	///< Number of total samples
			int children_level;	   	///< Level of children
			Node* child;	       	///< Pointer to child in the path (helps as to go down in the tree).
			Node* father;	       	///< Pointer to the father
			void* object;	       	///< the object stored
			bool descendant_flag;   ///< A new descendant has arrived in the family
			bool active_flag;	  	///< Indicates if the particular node is active or not (sampling procedure)
			int active_index;  		///< Indicates the index on the active nodes set.
			bool basis_flag;	   	///< Indicates if the particular node is an active basis or not
			int basis_index;
			real Q;                 ///< Value function
			
			///Bayesian Multivariate Model Prediction.
			real weight;    ///< backoff weight
			real weight_s;  ///< sampling weigth
			real posterior; ///< posterior probability
			BayesianMultivariateRegression* StatePrediction;
			BayesianMultivariateRegression* RewardPrediction;
			
			/// Constructor needs a point and a level
			Node(const CoverTree& tree_,const Vector& point_, const int level_, CoverTree::Node* const father_, void* object_);
			
			/// Destructor
			~Node();
			
			/// Metric
			const real distanceTo (const Vector& x) const
			{
				return tree.metric(x, point);
			}
			
			void UpdateStatistics(const Vector& state, const Vector& next_state, const real& reward, const bool& absorb, CoverTree::Node* const child_) 
			{
				Insert(state, next_state, reward, absorb);
				child = child_;
				samples++;
			}
			void Insert(const Vector& state, const Vector& next_state, const real& reward, const bool& absorb = false)
			{
				StatePrediction->AddElement(next_state, state);
				if(tree.RewardPred) {
					Vector r(reward);
					RewardPrediction->AddElement(r, state);
				}
			}
			real Update(const Vector& state, const Vector& next_state, real total_probability) {
				real p;
				real new_total_probability = 0;
				if(samples > 1) {
//					if(father == NULL)
//						printf("################################################root\n");
//				printf("Point\n");
//				state.print(stdout);
					p = StatePrediction->Posterior(state, next_state);
//				point.print(stdout);
//				printf("Childrens\n");
//				for(int i = 0; i<(int)children.size(); ++i) {
//					children[i]->point.printf(stdout);printf(" Samples %d",children[i]->samples);
//					printf("\n");
//				}
					real temp = weight*p;
					new_total_probability = temp + (1-weight)*total_probability;

//					real thres = 1.0 - 1e-6;
					real weight_new = temp / new_total_probability;
////
//					if(weight_new < 1e-6) {
//						weight_new = 1e-6;
//					}
//					if(weight_new > thres && father!=NULL) 
//						weight = thres;
//					else {
//						weight = weight_new;
//					}
					weight = weight_new;
	//				if(father == NULL)
//						weight = 1;
//					else 
//						weight = weight_new;
					
//					printf("probability = %f Samples = %d, weight = %f, level = %d, depth = %d\n",p,samples,weight,level,depth);

				}
				return new_total_probability;
			}
			
			/// Insert a new point at the given level, as a child of this node
			void Insert(const Vector& new_point, const int level, void* obj = NULL);
			void Insert(const Vector& new_point, const Vector& phi, const Vector& next_state, const real& reward, const bool& absorb, const int level, void* obj = NULL);
			
			/// Find nearest neighbour of the node
			std::pair<const CoverTree::Node*, real> NearestNeighbour(const Vector& query, const real distance) const;
			std::pair<CoverTree::Node*, real> FindNearestNeighbour(const Vector& query, const real distance);
			
			/// The number of children
			const int Size() const
			{
				return children.size();
			}
			/// A lower bound on the level of children
			const int ChildrenLevel() const
			{
				return children_level;
			}		
			///Find the number of descendants
			const int Descendants() const
			{
				return descendants;
			}		
			/// Returns the nodes' depth in the tree
			const int Depth() const
			{
				return depth;
			}
			/// Returns the number of observed points
			const int NumObs() const
			{
				return samples;
			}
			/// Change the sampling approach of node's model
			void SamplingModel(bool Thompson);
			/// Display the tree in textual format
			void Show() const;	
			/// Displaye the sampling tree in textual format
			void ShowSampling() const;
			void ShowBasis() const;
			/// Predict the next state given previous state.
			const Vector GenerateS(const Vector& state) const {
				return StatePrediction->generate(state);
			}
			/// Predict the next reward given previous state.
			const real GenerateR(const Vector& state) const {
				Vector r = RewardPrediction->generate(state);
				return r[0];
			}
			/// Returns the nodes weight
			const real GetWeight() const {
				return weight;
			}
			/// Returns the node's index
			const int GetIndex() const {
				return index;
			}
			const int GetSamplingIndex() const {
				return active_index;
			}
			const int GetBasisIndex() const {
				return basis_index;
			}
			const bool GetActive() const {
				return active_flag;
			}
			const bool GetActiveBasis() const {
				return basis_flag;
			}
			const real GetValueFunction() const {
				return Q;
			}
			const void SetValueFunction(const real& q_) {
				Q = q_;
			}
			/// Display the tree in .dot format
			void Show(FILE* fout) const;
			/// Displaye the sampling tree in .dot format
			void ShowSampling(FILE* fout) const;
			void ShowBasis(FILE* fout) const;
		};
		
		/** Cover set.
		 
		 This is primarily a simple encapsulation of a vector of Node*.
		 On the other hand, it also provides some set semantics:
		 a node value may not be duplicated within a set.
		 */
		struct CoverSet
		{
			std::vector<Node*> nodes;
			std::vector<real> distances;
			const int Size() const
			{
				return nodes.size();
			}
			void Insert(Node* node, real distance)
			{
				int N = Size();
				for (int i=0; i<N; ++i) {
					if (nodes[i] == node) {
						return;
					}
				}
				nodes.push_back(node);
				distances.push_back(distance);
			}
			const int NChildren(const int i) const
			{
				return nodes[i]->Size();
			}
			void Show() const 
			{
				for (int i=0; i<Size(); ++i) {
					printf ("%d : ", i);
					nodes[i]->Show();
				}
			}
		};
		
		const real	metric(const CoverSet& Q, const Vector& p) const;
		void		UpdateStatistics(const Vector& input, const Vector& output);
		const void	SamplingNode(Node* n);
		const void	SamplingNode(Node* n, const Vector& R);
		const void	SamplingTree();
		const void	SamplingTree(const Vector& R);
		
		Vector        BasisCreation(const Vector& state) const;
		const std::vector< std::pair<int,real> > ExternalBasisCreation(const Vector& state) const;
		
		Node*		Insert(const Vector& new_point, const CoverSet& Q_i, const int level, void* obj);
		Node*		Insert(const Vector& new_point, const Vector& phi, const Vector& next_state, const real& reward, const bool& absorb, const CoverSet& Q_i, const int level, void* obj);
		Node*		Insert(const Vector& new_point, void* obj = NULL);
		Node*		Insert(const Vector& new_point, const Vector& next_state, const real& reward, const bool& absorb = false, void* obj = NULL);
		
		const Node*  	NearestNeighbour(const Vector& query_point) const;
		const Node*  	SelectedNearestNeighbour(const Vector& query_point) const;
		const Vector	GenerateState(const Vector& query_point) const;
		const real   	GenerateReward(const Vector& query_point) const;
		const real    GetValueFunction(const Vector& query_point) const;
		const void    SetValueFunction(const Vector& query_point, const real& q);
		const int     GetNumNodes() const;
		const void   	SetNumNodes(int num);
		const int    	GetNumSamplingNodes() const;
		const void   	SetNumSamplingNodes(int num);
		const int    	GetNumBasisNodes() const;
		const void   	SetNumBasisNodes(int num); 
		const void   	Reset();
		const real   	GetEntranceThreshold(int depth) const;
		const real   	GetBasisThreshold() const;
		bool Check() const;
		void SamplingModel(bool Thompson);
		void Show() const;
		void ShowSampling() const;
		void ShowBasis() const;
		CoverTree(real c, real a_ = 0.1, real N0_ = 0.1, RBFBasisSet* RBFs_ = NULL, bool Sampling = true, bool f = false);
		~CoverTree();
	protected:
		bool	Check(const CoverSet& parents, const int level) const;
		real	Separation(const CoverSet& Q) const;
		Node*	FindNearestNeighbour(const Vector& query_point) const;
		int  	tree_level;
		int  	num_nodes;
		real	log_c;
		int  	thres;
		int  	num_sampling_nodes;
		int  	num_basis_nodes;
		int  	total_samples;		    
		///Regression parameters.
		real	a;
		real	N0;
		std::vector<Node*> Basis;
		////////////////////////
		RBFBasisSet* RBFs;
		bool RewardPred;		//Reward prediction
		bool ThompsonSampling;  //Thompson sampling
		Node* root;
	};


#endif
