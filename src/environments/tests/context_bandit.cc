#include "InventoryManagement.h"
#include "ValueIteration.h"
#include "Matrix.h"
#include "MatrixNorm.h"
#include "Gridworld.h"
#include "ContextBandit.h"
#include "QLearning.h"
#include "MersenneTwister.h"

#include <vector>

int main(void)
{
    real gamma = 0.9;
    int n_states = 8;
    int n_actions = 2;
    MersenneTwisterRNG rng; 
    ContextBandit context_bandit(n_states, n_actions, &rng);

    const DiscreteMDP* mdp = context_bandit.getMDP();
    ValueIteration value_iteration(mdp, gamma);
    value_iteration.ComputeStateActionValues(1e-6);
    Matrix Q = value_iteration.getValues();

    real lambda = 0.5;
    real alpha = 0.01;
    real epsilon = 0.1;
    VFExplorationPolicy* exploration_policy = 
        new EpsilonGreedy(n_actions, epsilon);
    QLearning q_learning(n_states, n_actions, gamma, lambda, alpha,
                         exploration_policy);
    int n_steps = 1000;
    for (int t=0; t<n_steps; ++t) {
        int state = context_bandit.getState();
        real reward = context_bandit.getReward();
        int action = q_learning.Act(reward, state);
        context_bandit.Act(action); // fail
    }

    DiscreteMDP* mdp_copy = new DiscreteMDP(*mdp);

    delete mdp_copy;
    delete mdp;

    
    delete exploration_policy;
    
    return 0;
}
