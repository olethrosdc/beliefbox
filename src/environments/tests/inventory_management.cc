#include "InventoryManagement.h"
#include "ValueIteration.h"
#include "MultiMDPValueIteration.h"
#include "Matrix.h"
#include "MatrixNorm.h"
#include "Gridworld.h"
#include "DiscreteMDP.h"
#include <vector>

Matrix GetQValues(const DiscreteMDP* mdp, real gamma)
{
    ValueIteration value_iteration(mdp, gamma);
    value_iteration.ComputeStateActionValues(10e-9);
    
    Matrix Q(mdp->getNStates(), mdp->getNActions());
    for (int s=0; s<mdp->getNStates(); ++s) {
        for (int a=0; a<mdp->getNActions(); ++a) {
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
    
    Matrix Q(M[0]->getNStates(), M[0]->getNActions());
    for (int s=0; s<M[0]->getNStates(); ++s) {
        for (int a=0; a<M[0]->getNActions(); ++a) {
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

#if 0
    std::vector<const DiscreteMDP*> mdp_list;
    std::vector<Matrix> Q;
    mdp_list.push_back(inventory_management_A.getMDP());
    mdp_list.push_back(inventory_management_B.getMDP());
    Q.push_back(GetQValues(mdp_list[0], gamma));
    Q.push_back(GetQValues(mdp_list[1], gamma));

    for (real weight = 0; weight<=1; weight+=0.1) {
        Vector w(2);
        w(0) = weight;
        w(1) = 1 - weight;
        Matrix Qw = GetQValues(mdp_list, w, gamma);
        printf ("%f %f\n", 
                FrobeniusNorm(Qw - Q[0]),
                FrobeniusNorm(Qw - Q[1]));
    }
#endif

#if 1
    int T = 1000;
    inventory_management_A.Reset();

    int n_actions = inventory_management_A.getNActions();
    for (int t=0; t<T; ++t) {
        int state = inventory_management_A.getState();
        //int reward = inventory_management_A.getReward();
        int a = rand()%n_actions;
        inventory_management_A.Act(a);
        
    }
#endif
    return 0;
}
