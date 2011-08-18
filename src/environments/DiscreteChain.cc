#include "DiscreteChain.h"
DiscreteChain::DiscreteChain(int n, real slip_, real start_, real end_)
    : Environment<int, int> (n, 2),
	  slip(slip_),
	  start(start_),
	  end(end_)
{
    logmsg ("Dscrete chain, slip: %f, r_s: %f, r_e: %f\n", slip, start, end);
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
	int forward = action;
	if (urandom() < slip) {
		forward = 1 - forward;
	}
	action = forward;
	switch(action) {
	case 0:
		reward = start;
		break;
	case 1:
		if (state == (int) n_states - 1) {
			reward = end;
		} else {
			reward = 0.0;
		}
 		break;
	}



	if (forward) {
		state++;
		if (state > (int) n_states - 1) {
			state = n_states - 1;
        }
	} else {
		state = 0;
	}

    return true;
}

DiscreteMDP* DiscreteChain::getMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);

    for (uint s=0; s<n_states; ++s) {
        // Action 0
        if (s==0) {
            mdp->setTransitionProbability(s, 0, 0, 1.0);
        } else {
            mdp->setTransitionProbability(s, 0, 0, 1.0 - slip);
        }
        for (uint j=1; j<n_states;++j) {
			if (j == s) {
                mdp->setTransitionProbability(s, 0, j, slip);
			} else {
                mdp->setTransitionProbability(s, 0, j, 0.0);
            }
        }

        // Action 1
        uint s_n = s + 1;
        if (s_n > n_states - 1) {
            s_n = n_states - 1;
        }
        mdp->setTransitionProbability(s, 1, s_n, 1 - slip);
        for (uint j=0; j<n_states;++j) {
			if (j == 0) {
                mdp->setTransitionProbability(s, 1, j, slip);
			} else if (j != s_n) {
                mdp->setTransitionProbability(s, 1, j, 0.0);
            }
        }
        if (s < n_states - 1) {
            mdp->addFixedReward(s, 0, start);
            mdp->addFixedReward(s, 1, 0.0);
        } else {
            mdp->addFixedReward(s, 0, start);
            mdp->addFixedReward(s, 1, end);
        }
    }
    
    return mdp;
}

