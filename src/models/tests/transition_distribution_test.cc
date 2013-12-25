// -*- Mode: c++ -*-

#include "TransitionDistribution.h"
#include "real.h"
#include "Random.h"

void DisplayTransitions(const DiscreteTransitionDistribution& kernel);

int main ()
{
	int n_states = 4;
	int n_actions = 1;
	real C = 0.1;
	DiscreteTransitionDistribution kernel(n_states, n_actions);

	printf("Adding transitions randomly\n");
	for (int i=0; i<n_states; i++) {
		for (int a=0; a<n_actions; a++) {
			int j = 0;
			real remaining = 1.0;
			while(remaining > 0) {
				j++;
				if (urandom() < C) {
					real p = remaining*urandom();
					remaining -= p;
					kernel.SetTransition(i, a, j, p);
				}
			}
		}
	}
	

	DisplayTransitions(kernel);

	printf("Clearing transitions for action 0\n");
	for (int i=0; i<n_states; i++) {
		int a = 0;
		for (int j=0; j<n_states; j++) {
			kernel.SetTransition(i, a, j, 0);
		}
	}

	DisplayTransitions(kernel);
	
	return 0;
}

void DisplayTransitions(const DiscreteTransitionDistribution& kernel) 
{
	
	printf("Displaying transitions\n");
	int n_states = kernel.GetNStates();
	int n_actions = kernel.GetNActions();
	for (int i=0; i<n_states; i++) {
		for (int a=0; a<n_actions; a++) {
			for (int j=0; j<n_states; j++) {
				real p = kernel.GetTransition(i, a, j);
				if (p > 0) {
					printf("%d %d %d %f\n", i, a, j, p);
				}
			}
		}
	}
}
