#include "InventoryManagement.h"
#include "ValueIteration.h"
#include "MultiMDPValueIteration.h"
#include "Matrix.h"
#include "Gridworld.h"
#include "DiscreteMDP.h"
#include <vector>

Matrix GetQValues(const DiscreteMDP* mdp, real gamma)
{
    ValueIteration value_iteration(mdp, gamma);
    value_iteration.ComputeStateActionValues(10e-9);
    
    Matrix Q(mdp->GetNStates(), mdp->GetNActions());
    for (int s=0; s<mdp->GetNStates(); ++s) {
        for (int a=0; a<mdp->GetNActions(); ++a) {
            Q(s,a) = value_iteration.getValue(s, a);
        }
    }
    return Q;
}

Matrix GetQValues(std::vector<const DiscreteMDP*>& M,
                  Vector& w,
                  real gamma)
{
    MultiMDPValueIteration value_iteration(w, M, gamma);
    value_iteration.ComputeStateActionValues(10e-9);
    
    Matrix Q(M[0]->GetNStates(), M[0]->GetNActions());
    for (int s=0; s<M[0]->GetNStates(); ++s) {
        for (int a=0; a<M[0]->GetNActions(); ++a) {
            Q(s,a) = value_iteration.getValue(s, a);
        }
    }
    return Q;
}


int main(void)
{
    real gamma = 0.9;
    InventoryManagement inventory_management_A(15, 15, 0.1, 1.1);
    InventoryManagement inventory_management_B(15, 15, 0.2, 1.1);

    std::vector<const DiscreteMDP*> mdp_list;
    std::vector<Matrix> Q;
    mdp_list.push_back(inventory_management_A.getMDP());
    mdp_list.push_back(inventory_management_B.getMDP());
    Q.push_back(GetQValues(mdp_list[0], gamma));
    Q.push_back(GetQValues(mdp_list[1], gamma));

    for (real weight = 0; weight<=1; weight+=0.01) {
        Vector w(2);
        w(0) = weight;
        w(1) = 1 - weight;
        Matrix Qw = GetQValues(mdp_list, w, gamma);
        printf ("%f %f\n", 
                (Qw - Q[0]).p_norm(2.0),
                (Qw - Q[1]).p_norm(2.0));
    }
    return 0;
}
