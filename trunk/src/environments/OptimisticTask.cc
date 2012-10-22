#include "OptimisticTask.h"

/// In the optimistic task, you get reward epsilon in the normal
/// state, reward 1 in the good state.
OptimisticTask::OptimisticTask(real epsilon_, real delta_) 
    : Environment<int, int> (3, 2),
      epsilon(epsilon_),
      delta(delta_)
{
    Reset();
}
    

void OptimisticTask::Reset()
{
    state = 0;
    reward = 0;
}

bool OptimisticTask::Act(int action)
{
    // get next state
    switch (state) {
    case 0:
        if (action == 0) {
            state = 0; 
        } else {
            state = 1;
        }
        break;
    case 1:
        if (urandom() < delta) {
            if (action == 0) {
                state = 0;
            } else {
                state = 2;
            }
        }
        break;
    case 2:
        if (action == 1) {
            state = 2;
        } else {
            state = 1;
        }
        break;
    }

    // get next reward
    switch(state) {
    case 0:
    case 2:
        reward = epsilon; 
        break;
    case 1:
        reward = 0;
        break;
    }

    return true;
}

DiscreteMDP* OptimisticTask::getMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);

    // get next state
    mdp->setTransitionProbability(0, 0, 0, 1);
    mdp->setTransitionProbability(0, 1, 1, 1);
    mdp->setTransitionProbability(2, 0, 1, 1);
    mdp->setTransitionProbability(2, 1, 2, 1);

    mdp->setTransitionProbability(1, 0, 0, delta);
    mdp->setTransitionProbability(1, 0, 1, 1 - delta);

    mdp->setTransitionProbability(1, 1, 2, delta);
    mdp->setTransitionProbability(1, 1, 1, 1 - delta);

    mdp->addFixedReward(0, 0, epsilon);
    mdp->addFixedReward(0, 1, epsilon);

    mdp->addFixedReward(1, 0, 0.0);
    mdp->addFixedReward(1, 1, 0.0);

    mdp->addFixedReward(2, 0, epsilon);
    mdp->addFixedReward(2, 1, epsilon);
    
    return mdp;
}

