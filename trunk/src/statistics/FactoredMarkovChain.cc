/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "FactoredMarkovChain.h"
#include "Random.h"

FactoredMarkovChain::FactoredMarkovChain(int n_actions_,
                                         int n_obs_,
                                         int mem_size_) 
    : n_actions(n_actions_),
      n_obs(n_obs_),
      mem_size(mem_size_),
      transitions((int) pow((double) (n_actions*n_obs), (double) mem_size), n_states),
      act_history(mem_size),
      obs_history(mem_size),
      history(mem_size),
      threshold(0.5)
{
    for (int i=0; i<mem_size; ++i) {
        act_history[i] = 0;
        obs_history[i] = 0;
        history[i] = 0;
    }
}

FactoredMarkovChain::~FactoredMarkovChain()
{
}

FactoredContext FactoredMarkovChain::CalculateStateID()
{
    MCState id = 0;
	MCState n = 1;
	for (MCState i=1; i<=mem_size; i++, n*=n_states) {
		id += history[i-1]*n; 
	}
	return id;
}
}

/** Train the chain with a new observation

   @return Probability of having observed the next state.
*/
real FactoredMarkovChain::Observe(int act, int prd)
{
    PushAction(act);
	curr_state = CalculateStateID();
	transitions.observe(curr_state, state);
    PushState (state);
    return Pr;
}


/** Probability of an observation given a particular action.

   \return The probability of observing prd given act.
*/
real FactoredMarkovChain::NextStateProbability (int act, int prd)
{
	assert((state>=0)&&(state<n_states));

	curr_state = CalculateStateID ();
	assert (curr_state>=0 && curr_state<tot_states);
	return getProbability(curr_state, act, prd);
}


void FactoredMarkovChain::getNextStateProbabilities(int act, std::vector<real>& p)
{
    curr_state = CalculateStateID();
    return getProbabilities(curr_state, p);
}



/// Get the number of transitions from \c ctx to \c prd
real FactoredMarkovChain::getTransition (FactoredContext ctx, int prd)
{
	assert((ctx>=0)&&(ctx<n_states));
	assert((prd>=0)&&(prd<n_states));
	return transitions.get_weight(ctx, prd);
}



/// Get the transition probability from \c ctx to \c prd
///
/// Takes into account the threshold.
real FactoredMarkovChain::getProbability(FactoredContext ctx, int prd)
{
	assert((ctx>=0)&&(ctx<tot_states));
	assert((prd>=0)&&(prd<n_states));
	real sum = 0.0;
	int N = transitions.nof_destinations();
	
	for (int i=0; i<N; ++i) {
		sum += transitions.get_weight(ctx, i);
	}
	
	return (transitions.get_weight(ctx, prd) + threshold) / (sum + threshold * ((real) N));
}

/// Get the transition probabilities
///
/// Takes into account the threshold.
void FactoredMarkovChain::getProbabilities(FactoredContext ctx, std::vector<real>& p)
{
	assert((ctx>=0)&&(ctx<tot_states));
	assert((int) p.size()== n_states);
	real sum = 0.0;
	int N = transitions.nof_destinations();
	
	for (int i=0; i<N; ++i) {
        p[i] = threshold + transitions.get_weight(ctx, i);
		sum += p[i];
	}

    real invsum = 1.0 / sum;
	for (int i=0; i<N; ++i) {
        p[i] *= invsum;
	}
}


/// Get the density of the transition probabilities \c p from \c ctx
real FactoredMarkovChain::pdf(FactoredContext ctx, Vector p)
{
	assert((ctx>=0)&&(ctx<tot_states));
	Swarning("Not implemented\n");
	return 0.0;
}




/**
   \brief Set the threshold specifying the Dirichlet prior.

   \arg \c threshold A threshold, added to all statistics. This way
   the probability of each transition will be \f$P(i|j)=\frac{N_{i,j}
   + \theta}{\sum_k N_{k,j} + \theta}\f$, there \f$\theta\f$ is the
   threshold.

   See also MarkovChainSoftmaxNormalize()
*/
void FactoredMarkovChain::setThreshold (real threshold)
{
	this->threshold = threshold;
}



/**
   \brief Reset the chain

   Sets the state (and history) to 0. Call before training for
   distinct sequences and before generating a new
   sequence. (Otherwise training and generation will depend on past
   sequences).

*/
void FactoredMarkovChain::Reset ()
{
	int i;
	for (i=0; i<mem_size; i++) {
		memory[i] = 0;
	}
	curr_state = 0;
}



/**
   \brief Generate values from the chain.

   Generates a new observation based on the current state history,
   according to the transition tables. Note that you might generate as
   many new states as you wish. The state history must be explicitly
   updated by calling MarkovChainPushState() with the generated value
   (or one of the generated values, if you have called this multiple
   times, or if you are selecting a generated value from a number of
   different markov chains).

   \arg \c chain: A pointer to the chain

   \return The ID of the next state, as generated by the MC. Returns -1 if
   nothing could be generated.

*/
int FactoredMarkovChain::GenerateStatic ()
{
	real tot = 0.0f;
        //int curr = CalculateStateID ();
	real sel = urandom();
	for (int j=0; j<n_states; j++) {
		real P = 0.0;//Pr[curr + j*tot_states];
		tot += P;
		if (sel<=tot && P>0.0f) {
			return j;
		}
	}

	return -1;
}

/** Add a new set of observations.

   Adds a new state to the history FIFO, pushing out the oldest
   member of the history. Note that the implementation could be
   faster, but the history is pretty small. 

   @param state A new state

   @return It returns the ID of the popped state.
*/
int  FactoredMarkovChain::PushState (int state) {
	int i;
	int popped_state;

	if (mem_size==0) { 
		return 0;
	}
	popped_state = memory[mem_size - 1];
	for (i = mem_size - 1; i>0; i--) {
		memory[i] = memory[i-1];
	}
	memory[0] = state;
    
	//logmsg("Pushing %d, popping %d\n", state, popped_state);
	return popped_state;
}

/** Can be useful for debugging.
*/
int FactoredMarkovChain::ShowTransitions () {
	printf ("\nState transition dump %ld\n", tot_states);
	for (int i=0; i<tot_states; i++) {
		printf ("Transition %d : ", i);
        std::vector<real> p(n_states);
        getProbabilities(i, p);
		for (int j=0; j<n_states; j++) {
            printf (" %f", p[j]);
		}
        printf("\n");
	}
	return 0;
}

