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

#include "KNNModel.h"
#include "BasisSet.h"
#include "Distribution.h"
#include "Random.h"

KNNModel::KNNModel(int n_actions_, int n_dim_, real gamma_, bool optimistic, real optimism_, real r_max_)
    : n_actions(n_actions_), n_dim(n_dim_), kd_tree(n_actions_), gamma(gamma_),
      optimistic_values(optimistic), optimism(optimism_), r_max(r_max_),
      max_samples(-1)
{
    for (int i=0; i<n_actions; ++i) {
        kd_tree[i] = new KDTree<TrajectorySample> (n_dim);
    }
}

KNNModel::~KNNModel()
{
    for (int i=0; i< n_actions; ++i) {
        delete kd_tree[i];
    }
}

void KNNModel::AddSample(TrajectorySample x)
{
    if (max_samples < 0 || samples.size() < (uint) max_samples) {
        samples.push_back(x);
        kd_tree[x.a]->AddVectorObject(x.s, &samples.back());
        x.dV = 1.0;
        x.V = x.r;
    }
}



/** Predict the next reward and a state
    
    This is the simplest type of predictor.
    It just predicts the next vector mean, no distribution is used.
    
    \param alpha absolute-to-relative parameter. Set to 0 for absolute predictions, 1 for predictions relative to the current state.
    \param x starting state.
    \param action action
    \param reward reward reference (to return).
    \param y end state reference (to return).
    \param K numbers of neighbours to consider.
    \param b RBF width to use.
*/
void KNNModel::GetExpectedTransition(real alpha, Vector& x, int action, real& reward, Vector& y, int K, real b)
{

    RBF rbf(x, b);

    //basis.Evaluate(x);
    assert(n_dim == y.Size());
    reward = 0.0;
    for (int i=0; i<n_dim; ++i) {
        y[i] = 0;
    }

    OrderedFixedList<KDNode> node_list = kd_tree[action]->FindKNearestNeighbours(x, K);
    
    std::list<std::pair<real, KDNode*> >::iterator it;
    
    real sum = 0;
    Vector w(K);
    int i=0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it, ++i) {
        KDNode* node = it->second;
        TrajectorySample* sample = kd_tree[action]->getObject(node);
        w[i] =  rbf.Evaluate(sample->s);
        sum += w[i];
    }
    w /= sum;
    i = 0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it, ++i) {
        KDNode* node = it->second;
        TrajectorySample* sample = kd_tree[action]->getObject(node);
        y += (sample->s2 + (x - sample->s)*alpha)* w[i];
        reward += sample->r * w[i];
    }
}

/** Get the expected value for an action
    
    Add the probability that there is a posibility to go to an 
    unknown terminal state with highish value.
 */
real KNNModel::GetExpectedActionValue(Vector& x, int action, int K, real b)
{
    RBF rbf(x, b);

    OrderedFixedList<KDNode> node_list = kd_tree[action]->FindKNearestNeighbours(x, K);
    
    std::list<std::pair<real, KDNode*> >::iterator it;
    
    real sum = 0;
    real Q = 0.0;
    for (it = node_list.S.begin(); it != node_list.S.end(); ++it) {
        KDNode* node = it->second;
        TrajectorySample* sample = kd_tree[action]->getObject(node);
        real w =  rbf.Evaluate(sample->s);
        Q += sample->V * w;
        sum += w;
    }
    if (optimistic_values) {
        Q += optimism * r_max / (1.0 - gamma);
        sum += optimism;
    }
    Q /= sum;
    if (isnan(Q)) {
        Q = 0.0;
    }
    return Q;
}

// Return the expected value of a state according to our model
real KNNModel::GetExpectedValue(Vector& x, int K, real b)
{
    Vector Q(n_actions);
    for (int a=0; a<n_actions; a++) {
        Q[a] = GetExpectedActionValue(x, a, K, b);

    }
    return Max(&Q);
}

// Return the expected value of a state according to our model
int KNNModel::GetBestAction(Vector& x, int K, real b)
{

    Vector Q(n_actions);
    Vector p(n_actions);
    for (int a=0; a<n_actions; a++) {
        Q[a] = GetExpectedActionValue(x, a, K, b);
    }
    //SoftMax(Q.Size(), &Q[0], &p[0], greedy);
    //return DiscreteDistribution::generate(p);
    return ArgMax(&Q);
}


/// Update the value of the current state
void KNNModel::UpdateValue(TrajectorySample& start_sample, real alpha, int K, real b)
{
    if (start_sample.terminal) {
        start_sample.V = start_sample.r;
        start_sample.dV = 0.0;
        return;
    }
    
    RBF rbf(start_sample.s, b);
    Vector Q(n_actions);

    for (int a=0; a<n_actions; ++a) {
        OrderedFixedList<KDNode> node_list = kd_tree[a]->FindKNearestNeighbours(start_sample.s, K);
        
        std::list<std::pair<real, KDNode*> >::iterator it;
    
        real sum = 0;
        Q[a] = 0.0;
        //printf("Action %d: ", a);
        for (it = node_list.S.begin(); it != node_list.S.end(); ++it) {
            KDNode* node = it->second;
            TrajectorySample* sample = kd_tree[a]->getObject(node);
            Vector y = sample->s2 + (start_sample.s - sample->s) * alpha;
            real w =  rbf.Evaluate(sample->s);
            real Qa_i = (sample->r + gamma*GetExpectedValue(y, K, b));
            if (isnan(Qa_i)) {
                Qa_i = 0.0;
            }
            //printf ("%f (%f) ", Qa_i, w);
            Q[a] += Qa_i*w;
            sum += w;
        }
        Q[a] /= sum;
        if (node_list.S.size() == 0) {
            Q[a] = 0.0;
        }
        //printf ("-> %f\n", Q[a]);
    }

    real new_V = Max(&Q);
    if (isnan(new_V)) {
        new_V = 0.0;
    }
    start_sample.dV = fabs(start_sample.V - new_V);
    start_sample.V = new_V;
    
}

/** Perform value iteration
    
    Iterate over all states and update their values
 */

void KNNModel::ValueIteration(real alpha, int K, real b)
{
    for (std::list<TrajectorySample>::iterator it = samples.begin();
         it != samples.end(); ++it) {
        if (it->dV > 0) {//10e-6) {
            UpdateValue(*it, alpha, K, b);
        }
    }
    
}


void KNNModel::Show()
{
    for (std::list<TrajectorySample>::iterator it = samples.begin();
         it != samples.end(); ++it) {
        for (int i=0; i<n_dim; ++i) {
            printf ("%f ", it->s[i]);
        }
        printf ("%d %f ", it->a, it->r);
        for (int i=0; i<n_dim; ++i) {
            printf ("%f ", it->s2[i]);
        }
        printf ("%f %d # SAMPLE\n", it->V, it->terminal);
    }
    
}
