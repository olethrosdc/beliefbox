#include "DiscreteChain.h"
DiscreteChain::DiscreteChain(int n, real slip_, real start_, real end_)
    : Environment<int, int> (n, 2),
	  slip(slip_),
	  start(start_),
	  end(end_)
{
  //logmsg ("Discrete chain, slip: %f, r_s: %f, r_e: %f\n", slip, start, end);
    mdp = getMDP();
    Reset();

}
    
DiscreteChain::~DiscreteChain()
{
    delete mdp;
}

void DiscreteChain::Reset()
{
    state = 0;
    reward = 0;
    mdp->Reset(0);
}

bool DiscreteChain::Act(const int& action)
{
    int action_taken = action;
	int forward = action;
	if (urandom() < slip) {
		forward = 1 - forward;
	}
	action_taken = forward;
	switch(action_taken) {
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
        uint s_n = s + 1;
        if (s_n > n_states - 1) {
            s_n = n_states - 1;
        }
        // Action 0
        if (s==0) {
            mdp->setTransitionProbability(s, 0, 0, 1.0);
        } else {
            mdp->setTransitionProbability(s, 0, 0, 1.0 - slip);
	    mdp->setTransitionProbability(s, 0, s_n, slip);
        }

#if 1
        for (uint j=1; j<n_states;++j) {
	  if (j != s_n && j != 0) {
                mdp->setTransitionProbability(s, 0, j, 0.0);
            }
        }
#endif

        // Action 1
        mdp->setTransitionProbability(s, 1, s_n, 1 - slip);
        for (uint j=0; j<n_states;++j) {
			if (j == 0) {
                mdp->setTransitionProbability(s, 1, j, slip);
			} else if (j != s_n) {
                mdp->setTransitionProbability(s, 1, j, 0.0);
            }
        }
        if (s < n_states - 1) {
            mdp->addFixedReward(s, 0, start * (1 - slip));
            mdp->addFixedReward(s, 1, start * slip);
        } else {
            mdp->addFixedReward(s, 0, start * (1 - slip) + slip * end);
            mdp->addFixedReward(s, 1, end * (1 - slip) + slip * start);
        }
    }
    
    return mdp;
}

