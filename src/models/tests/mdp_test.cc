// -*- Mode: c++ -*-

#ifdef MAKE_MAIN

#include "DiscreteMDPCounts.h"
#include "DiscreteChain.h"

int main (void)
{
	int n_states = 5;
	int n_actions = 4;
	DiscreteMDPCounts model(n_states, n_actions);
	DiscreteChain chain(5);
	DiscreteEnvironment* environment = &chain;
	for (int i=0; i<10; ++i) {
		int state = environment->getState();
		int action = rand()%n_actions;
		environment->Act(action);
		real reward = environment->getReward();
		int next_state = environment->getState();
		model.AddTransition(state, action, reward, next_state);
	}

	printf("Checking Mean MDP\n");
    const DiscreteMDP* mean_mdp = model.getMeanMDP();
	mean_mdp->Check();

	int n_tests = 10;
	for (int i=0; i<n_tests; ++i) {
		printf("Checking MDP %d/%d\n", i, n_tests);
		DiscreteMDP* sample_mdp = model.generate();
		sample_mdp->Check();
	}
	
	return 0;
}

#endif
