#include "KDTree.h"
#include "Vector.h"
KDTree::KDTree(int n) : n_dimensions(n), box_sup(n), box_inf(n), root(NULL)
{
    for(int i=0; i<n; ++i) {
        box_sup[i] = RAND_MAX;
        box_inf[i] = -RAND_MAX;
    }
}
KDTree::~KDTree()
{
    int N = node_list.size();
    for (int i=0; i<N; ++i) {
        delete node_list[i];
    }
}
void KDTree::AddVector(Vector& x)
{
    if (!root) {
        root = new KDNode(x, 0, box_inf, box_sup);
        node_list.push_back(root);
        return;
    }
    KDNode* node = root->AddVector(x, box_inf, box_sup);
    if (node) {
        node_list.push_back(node);
    } else{
        fprintf(stderr, "Strange: node not added\n");
    }
}


void KDTree::Show()
{
    int N = node_list.size();
    for (int i=0; i<N; ++i) {
        KDNode* node = node_list[i];
        printf ("N[%d] = {", i);
        for (int j=0; j<n_dimensions; ++j) {
            printf ("%.2f ", node->c[j]);
        }
        printf(" | %d}, [ (", node->a);
        
        for (int j=0; j<n_dimensions; ++j) {
            printf ("%.2f ", node->box_inf[j]);
        }
        printf("), (");
        for (int j=0; j<n_dimensions; ++j) {
            printf ("%.2f ", node->box_sup[j]);
        }
        printf(")]\n");
    }
}
    
KDNode* KDTree::FindNearestNeighbourLinear(Vector& x)
{
    int N = node_list.size();
    real min_dist = L1Norm(&x, &node_list[0]->c);
    KDNode* arg_min = node_list[0];
    for (int i=1; i<N; ++i) {
        real dist = L1Norm(&x, &node_list[i]->c);
        if (dist < min_dist) {
            min_dist = dist;
            arg_min = node_list[i];
        }
    }
    return arg_min;
}

KDNode* KDTree::FindNearestNeighbour(Vector& x)
{
    if (!root) {
        // tree is empty
        return NULL; 
    }

    KDNode* node = root->NearestNeighbour(x);
    return node;
}


KDNode* KDNode::AddVector(Vector& x, Vector& inf, Vector& sup)
{
    box_inf = inf;
    box_sup = sup;

    
    // if we already have a child, just go down
    if (x[a] <= c[a]) {
        sup[a] = c[a];
        if (lower) {
            return lower->AddVector(x, inf, sup);
        } else {
            Vector diff = sup - inf;
            lower = new KDNode(x, ArgMax(&diff), inf, sup);
            return lower;
        }
    } else {
        inf[a] = c[a];
        if (upper) {
            return upper->AddVector(x, inf, sup);
        } else {
            Vector diff = sup-inf;
            upper = new KDNode(x, ArgMax(&diff), inf, sup);
            return upper;
        }
    }
    return NULL;
}

KDNode* KDNode::NearestNeighbour(Vector& x, real dist)
{
    // if we already have a child, just go down
    if (x[a] <= c[a]) {
        if (lower) {
            return lower->NearestNeighbour(x);
        } else {
            return this;
        }
    } else {
        if (upper) {
            return upper->NearestNeighbour(x);
        } else {
            return this;
        }
    }
    return NULL;
}
