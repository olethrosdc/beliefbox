
#include "RiverSwim.h"
RiverSwim::RiverSwim(int n, real r_start_, real r_end_)
    : Environment<int, int> (n + 1, 2),
	  r_start(r_start_),
	  r_end(r_end_),
	  p_right(0.3),
	  p_left(0.1),
	  p_stuck(1.0 - (p_right + p_left))
{
	model = getMDP();
	model->Check();
	Reset();
}
    
RiverSwim::~RiverSwim()
{
	delete model;
}

/// Start in state 1 or 2 (2 or 3 here) with equal probability
void RiverSwim::Reset()
{
	if (urandom() < 0.5) {
		state = 2;
	} else {
		state = 3;
	}
    reward = 0;
	model->setState(state);
}

bool RiverSwim::Act(const int& action)
{
#if 0
	for (uint i=0; i<n_states; ++i) {
		printf ("%f ", model->getTransitionProbability(state, action, i));
	}
	printf("= P(s' | s = %d, a = %d)\n", state, action);
#endif
	reward = model->Act(action);
	state = model->getState();
    return true;
}

/** getMDP - return a pointer to a newly-created MPD structure

	This is a modification of the actual RiverSwim so that it
	corresponds properly to the state-action-state reward model it
	uses. Thus, two further states are added.
 */
DiscreteMDP* RiverSwim::getMDP() const
{
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);
	
	// -- intermediate states --
    for (uint s=1; s<n_states - 2; ++s) {
		// state-basd rewards
		mdp->addFixedReward(s, 0, 0.0);
		mdp->addFixedReward(s, 1, 0.0);

        // Action 0 always goes left
        mdp->setTransitionProbability(s, 0, s - 1, 1.0);
		
        // Action 1 sometimes goes left, or gets stuck
        mdp->setTransitionProbability(s, 1, s + 1, p_right);
        mdp->setTransitionProbability(s, 1, s, p_stuck);
        mdp->setTransitionProbability(s, 1, s - 1, p_left);
    }

	// -- beginning state --
	// state-based rewards
	mdp->addFixedReward(0, 0, r_start);
	mdp->addFixedReward(0, 1, 0.0);
	mdp->setTransitionProbability(0, 0, 0, 1.0);
	mdp->setTransitionProbability(0, 1, 0, 1.0 - p_right);
	mdp->setTransitionProbability(0, 1, 1, p_right);

	// Action 0 always goes left
	mdp->setTransitionProbability(n_states - 2, 0, n_states - 3, 1.0);
	mdp->setTransitionProbability(n_states - 1, 0, n_states - 3, 1.0);

	// Action 1 might go right!
	mdp->setTransitionProbability(n_states - 2, 1, n_states - 1, p_right);
	mdp->setTransitionProbability(n_states - 2, 1, n_states - 3, 1 - p_right);
	mdp->setTransitionProbability(n_states - 1, 1, n_states - 1, p_right);
	mdp->setTransitionProbability(n_states - 1, 1, n_states - 3, 1 - p_right);

	// -- final state --
	// state-based rewards
	mdp->addFixedReward(n_states - 2, 0, 0.0);
	mdp->addFixedReward(n_states - 2, 1, 0.0);
	mdp->addFixedReward(n_states - 1, 0, r_end);
	mdp->addFixedReward(n_states - 1, 1, r_end);

    return mdp;
}

