#include "DiscreteChain.h"
DiscreteChain::DiscreteChain(int n) 
    : Environment<int, int> (n, 2)
{
    
}
    
DiscreteChain::~DiscreteChain()
{
}

void DiscreteChain::Reset()
{
    state = 0;
    reward = 0;
}

bool DiscreteChain::Act(int action)
{
    if (state == (int) n_states - 1) {
        reward = 1;
    } else {
        reward = 0;
    }

    if (action == 0) {
        state = 0;
    }

    if (action == 1) {
        state++;
        if (state > (int) n_states - 1) {
            state = n_states - 1;
        }
    }
    return true;
}

DiscreteMDP* DiscreteChain::getMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);

    for (int s=0; s<n_states; ++s) {
        // Action 0
        mdp->setTransitionProbability(s, 0, 0, 1.0);
        for (int j=0; j<n_states;++j) {
            if (j != 0) {
                mdp->setTransitionProbability(s, 0, j, 0.0);
            }
        }

        // Action 1
        int s_n = s + 1;
        if (s_n > n_states - 1) {
            s_n = n_states - 1;
        }
        mdp->setTransitionProbability(s, 1, s_n, 1);
        for (int j=0; j<n_states;++j) {
            if (j != s_n) {
                mdp->setTransitionProbability(s, 1, j, 0.0);
            }
        }
        if (s < n_states - 1) {
            mdp->addFixedReward(s, 0, 0.0);
            mdp->addFixedReward(s, 1, 0.0);
        } else {
            mdp->addFixedReward(s, 0, 1.0);
            mdp->addFixedReward(s, 1, 1.0);
        }
    }
    
    return mdp;
}

