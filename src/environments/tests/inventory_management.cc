#include "InventoryManagement.h"
#include "ValueIteration.h"
#include "Matrix.h"

int main(void)
{
    real gamma = 0.5;
    InventoryManagement inventory_management(5, 5, 0.1, 1.1);
    const DiscreteMDP* mdp = inventory_management.getMDP();
    
    mdp->ShowModel();
    ValueIteration value_iteration(mdp, gamma);
    value_iteration.ComputeStateActionValues(10e-6,0);
    value_iteration.ComputeStateValues(10e-6,0);

    Matrix Q(mdp->GetNStates(), mdp->GetNActions());
    for (int s=0; s<mdp->GetNStates(); ++s) {
        for (int a=0; a<mdp->GetNActions(); ++a) {
            Q(s,a) = value_iteration.getValue(s, a);
        }
        printf("V[%d] = %f\n", s, value_iteration.getValue(s));
    }
    Q.print(stdout);

    return 0;
}
