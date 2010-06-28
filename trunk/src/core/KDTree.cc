#include "KDTree.h"
#include "Vector.h"
/// Create a tree
void_KDTree::void_KDTree(int n) : n_dimensions(n), box_sup(n), box_inf(n), root(NULL)
{
    for(int i=0; i<n; ++i) {
        box_sup[i] = RAND_MAX;
        box_inf[i] = -RAND_MAX;
    }
}
/// Destroy the tree
void_KDTree::~void_KDTree()
{
    int N = node_list.size();
    for (int i=0; i<N; ++i) {
        delete node_list[i];
    }
}

/// Add a vector and associated object, creating a node on the fly.
void void_KDTree::AddVector(Vector& x, void* object)
{
    if (!root) {
        root = new KDNode(x, 0, box_inf, box_sup, object);
        node_list.push_back(root);
        return;
    }
    KDNode* node = root->AddVector(x, box_inf, box_sup, object);
    if (node) {
        node_list.push_back(node);
    } else{
        fprintf(stderr, "Strange: node not added\n");
    }
}

/// Display the tree
void void_KDTree::Show()
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
    
/// Find the nearest neighbour to x by linear search
KDNode* void_KDTree::FindNearestNeighbourLinear(Vector& x)
{
    int N = node_list.size();
    real min_dist = EuclideanNorm(&x, &node_list[0]->c);
    KDNode* arg_min = node_list[0];
    for (int i=1; i<N; ++i) {
        real dist = EuclideanNorm(&x, &node_list[i]->c);
        if (dist < min_dist) {
            min_dist = dist;
            arg_min = node_list[i];
        }
    }
    return arg_min;
}

/// Find the nearest neighour ro x
KDNode* void_KDTree::FindNearestNeighbour(Vector& x)
{
    if (!root) {
        // tree is empty
        return NULL; 
    }
    real dist = RAND_MAX;
    KDNode* node = NULL;
    root->NearestNeighbour(x, node, dist);

    return node;
}


/// Find the first K nearest neigbours to X by linear search through the list.
OrderedFixedList<KDNode> void_KDTree::FindKNearestNeighboursLinear(Vector& x, int K)
{
    int N = node_list.size();
    OrderedFixedList<KDNode> knn_list(K);
    for (int i=0; i<N; ++i) {
        real dist = EuclideanNorm(&x, &node_list[i]->c);
        knn_list.AddPerhaps(dist, node_list[i]);
    }
    return knn_list;
}

/// Find the first K nearest neighbours to x.
OrderedFixedList<KDNode> void_KDTree::FindKNearestNeighbours(Vector& x, int K)
{
    OrderedFixedList<KDNode> knn_list(K);
    if (!root) {
        return knn_list;
    }
    real dist = RAND_MAX;
    root->KNearestNeighbours(x, knn_list, dist);

    return knn_list;
}


/// Add a point to the corresponding (upper or lower) half, creating it if necessary.
KDNode* KDNode::AddVector(Vector& x, Vector& inf, Vector& sup, void* object)
{
    box_inf = inf;
    box_sup = sup;

    
    // if we already have a child, just go down
    if (x[a] <= c[a]) {
        sup[a] = c[a];
        if (lower) {
            return lower->AddVector(x, inf, sup, object);
        } else {
            Vector diff = sup - inf;
            lower = new KDNode(x, ArgMax(&diff), inf, sup, object);
            return lower;
        }
    } else {
        inf[a] = c[a];
        if (upper) {
            return upper->AddVector(x, inf, sup, object);
        } else {
            Vector diff = sup-inf;
            upper = new KDNode(x, ArgMax(&diff), inf, sup, object);
            return upper;
        }
    }
    return NULL;
}

/** Find nearest neighbour

    At each node \f$i\f$, corresponding to the subset \f$S_i\f$, find
    all children \f$j \in C(i)\f$ with distance bound \f$D_L(x, S_j) < d\f$.
    Look at the child with the smallest bound first.
    
    The distance bound is easy to calculate. If we are above \f$c_a\f$
    then the upper subset has bound 0, while the lower partition has
    bound \f$c_a - x_a\f$, since in the best case, the closest point
    \f$y\f$ in that subset will have \f$y_k = x_k\f$ for all 
    \f$k \neq a\f$ and \f$y_a = c_a\f$.
 */
void KDNode::NearestNeighbour(Vector& x, KDNode*& nearest, real& dist)
{
    real c_dist = EuclideanNorm(&x, &c);
    if (c_dist < dist) {
        nearest = this;
        dist = c_dist;
    }
    real delta = x[a] - c[a];
    
    // check the set with the lowest bound first
    KDNode* first;
    KDNode* second;

    if (delta <= 0) {
        first = lower;
        second = upper;
    } else {
        first = upper;
        second = lower;
    }
        
    // the first always has potential to reduce dist
    if (first) {
        first->NearestNeighbour(x, nearest, dist);
    }
    // only check second if dist has a chance to be reduced
    if ((fabs(delta) < dist) && second) {
        second->NearestNeighbour(x, nearest, dist);
    }

}

/** Find K nearest neighbour s

    This is a simple generalization NearestNeighbour().
    Instead of having a fixed distance bound, we, in addition,
    have a dynamic bound based on the K-th neighbour found so far.
 */
void KDNode::KNearestNeighbours(Vector& x,
                               OrderedFixedList<KDNode>& knn_list,
                               real& dist)
{
    real c_dist = EuclideanNorm(&x, &c);
    if (knn_list.AddPerhaps(c_dist, this)) {
        dist = std::min(dist, knn_list.UpperBound());
    }
    real delta = x[a] - c[a];
    
    // check the set with the lowest bound first
    KDNode* first;
    KDNode* second;

    if (delta <= 0) {
        first = lower;
        second = upper;
    } else {
        first = upper;
        second = lower;
    }
        
    // the first always has potential to reduce dist
    if (first) {
        first->KNearestNeighbours(x, knn_list, dist);
    }
    // only check second if dist has a chance to be reduced
    if ((fabs(delta) < dist) && second) {
        second->KNearestNeighbours(x, knn_list, dist);
    }

}
