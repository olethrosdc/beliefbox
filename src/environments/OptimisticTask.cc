#include "OptimisticTask.h"

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
    if (action == 0) {
        state = 0;
    } else { // action 1
        if (state < 2) {
            if (urandom() < delta) {
                state = 2;
            } else {
                state = 1;
            }
        }  else {
            state = 2;
        }
        
    }

    switch(state) {
    case 0: reward = epsilon; break;
    case 1: reward = 0; break;
    case 2: reward = 1; break;
    }

    return true;
}

DiscreteMDP* OptimisticTask::getMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);

    int good_state = 2;
    int bad_state = 1;
    int default_state = 0;
    // action 0 takes you to state 0 always
    for (int s=0; s<3; ++s) {
        mdp->setTransitionProbability(s, 0, default_state, 1);
    }

    // action 1 takes you to state 2 w.p. delta, and to state 1 as long
    // as you are in one of the first two states.
    for (int s=0; s<2; ++s) {
        mdp->setTransitionProbability(s, 1, good_state, delta);
        mdp->setTransitionProbability(s, 1, bad_state, 1.0 - delta);
    }
    // you can stay in the good state forever by taking action 1!
    mdp->setTransitionProbability(good_state, 1, good_state, 1.0);
    
    // default state
    mdp->addFixedReward(default_state, 0, epsilon);
    mdp->addFixedReward(default_state, 1, epsilon);

    // 'bad' state
    mdp->addFixedReward(bad_state, 0, 0.0);
    mdp->addFixedReward(bad_state, 1, 0.0);

    // 'good' state
    mdp->addFixedReward(good_state, 0, 1.0);
    mdp->addFixedReward(good_state, 1, 1.0);
    
    return mdp;
}

