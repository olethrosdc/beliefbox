
#include "DoubleLoop.h"
DoubleLoop::DoubleLoop(real r_left_, real r_right_)
    : Environment<int, int> (9, 2),
	  r_left(r_left_),
	  r_right(r_right_)
{
	model = getMDP();
	model->Check();
	Reset();
}
    
DoubleLoop::~DoubleLoop()
{
	delete model;
}

/// Start in state 1 or 2 (2 or 3 here) with equal probability
void DoubleLoop::Reset()
{
    state = 0;
    reward = 0;
	model->setState(state);
}

bool DoubleLoop::Act(const int& action)
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

	This is a modification of the actual DoubleLoop so that it
	corresponds properly to the state-action-state reward model it
	uses. Thus, two further states are added.
 */
DiscreteMDP* DoubleLoop::getMDP() const
{
   
    DiscreteMDP* mdp = new DiscreteMDP(n_states, n_actions);
    mdp->addFixedReward(0, 0, 0.0);
    mdp->addFixedReward(0, 1, 0.0);
	mdp->setTransitionProbability(0, 1, 5, 1.0);
	mdp->setTransitionProbability(0, 0, 1, 1.0);

    mdp->addFixedReward(1, 0, 0.0);
    mdp->addFixedReward(1, 1, 0.0);
    mdp->addFixedReward(2, 0, 0.0);
    mdp->addFixedReward(2, 1, 0.0);
    mdp->addFixedReward(3, 0, 0.0);
    mdp->addFixedReward(3, 1, 0.0);
    mdp->addFixedReward(4, 0, r_right);
    mdp->addFixedReward(4, 1, r_right);

    mdp->addFixedReward(5, 0, 0.0);
    mdp->addFixedReward(5, 1, 0.0);
    mdp->addFixedReward(6, 0, 0.0);
    mdp->addFixedReward(6, 1, 0.0);
    mdp->addFixedReward(7, 0, 0.0);
    mdp->addFixedReward(7, 1, 0.0);
    mdp->addFixedReward(8, 0, r_left);
    mdp->addFixedReward(8, 1, r_left);

	mdp->setTransitionProbability(1, 0, 2, 1.0);
	mdp->setTransitionProbability(2, 0, 3, 1.0);
	mdp->setTransitionProbability(3, 0, 4, 1.0);
	mdp->setTransitionProbability(4, 0, 0, 1.0);
	mdp->setTransitionProbability(5, 0, 0, 1.0);
	mdp->setTransitionProbability(6, 0, 0, 1.0);
	mdp->setTransitionProbability(7, 0, 0, 1.0);
	mdp->setTransitionProbability(8, 0, 0, 1.0);

	mdp->setTransitionProbability(5, 1, 6, 1.0);
	mdp->setTransitionProbability(6, 1, 7, 1.0);
	mdp->setTransitionProbability(7, 1, 8, 1.0);
	mdp->setTransitionProbability(8, 1, 0, 1.0);
	mdp->setTransitionProbability(1, 1, 2, 1.0);
	mdp->setTransitionProbability(2, 1, 3, 1.0);
	mdp->setTransitionProbability(3, 1, 4, 1.0);
	mdp->setTransitionProbability(4, 1, 0, 1.0);



    return mdp;
}

