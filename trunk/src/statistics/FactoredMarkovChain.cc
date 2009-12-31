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
    : T(0), 
      n_actions(n_actions_),
      n_obs(n_obs_),
      n_states(n_obs*n_actions),
      mem_size(mem_size_),
      transitions((int) pow((double) n_states, (double) mem_size), n_obs),
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

/** Calculate the current context.

    Calculates the current contexts up to the time the
    last action was taken.
 */
FactoredMarkovChain::Context FactoredMarkovChain::CalculateContext()
{
    Context context = 0;
	Context n = 1;
	for (Context i=1; i<=mem_size; i++, n*=n_states) {
		context += history[i-1]*n;
	}
	return context;
}

FactoredMarkovChain::Context FactoredMarkovChain::CalculateContext(int act)
{
    Context context = 0;
	Context n = 1;
	for (Context i=1; i<=mem_size; i++, n*=n_states) {
		context += history[i-1]*n;
	}
	return context;
}


/** Train the chain with a new observation
    
    A side-effect is that it updates the current_context.

    @return Probability of having observed the next state.
*/
real FactoredMarkovChain::Observe(int act, int prd)
{
    PushAction(act);
	current_context = CalculateContext();
    real Pr = getProbability(current_context, prd);
	transitions.observe(current_context, prd);
    PushObservation(prd);
    return Pr;
}



/** Probability of a next observation

   \return The probability of observing x
*/
real FactoredMarkovChain::ObservationProbability (int x)
{
	return getProbability(current_context, x);
}

/** Probability of an observation given a particular action.

   \return The probability of observing prd given act.
*/
real FactoredMarkovChain::ObservationProbability (int a, int x)
{
	assert((a>=0)&&(a<n_acts));
    
	return getProbability(getContext(a), x);
}

#if 0
void FactoredMarkovChain::getNextStateProbabilities(int act, std::vector<real>& p)
{
    curr_state = CalculateStateID();
    return getProbabilities(curr_state, p);
}
#endif


/// Get the number of transitions from \c context to \c prd
real FactoredMarkovChain::getTransition (Context context, int prd)
{
	assert((context>=0)&&(context<n_states));
	assert((prd>=0)&&(prd<n_states));
	return transitions.get_weight(context, prd);
}



/// Get the transition probability from \c context to \c prd
///
/// Takes into account the threshold.
real FactoredMarkovChain::getProbability(Context context, int prd)
{
	assert((context>=0)&&(context<tot_states));
	assert((prd>=0)&&(prd<n_states));
	real sum = 0.0;
	int N = transitions.nof_destinations();
	
	for (int i=0; i<N; ++i) {
		sum += transitions.get_weight(context, i);
	}
	
	return (transitions.get_weight(context, prd) + threshold) / (sum + threshold * ((real) N));
}

/// Get the transition probabilities
///
/// Takes into account the threshold.
void FactoredMarkovChain::getProbabilities(Context context, std::vector<real>& p)
{
	assert((context>=0)&&(context<tot_states));
	assert((int) p.size()== n_states);
	real sum = 0.0;
	int N = transitions.nof_destinations();
	
	for (int i=0; i<N; ++i) {
        p[i] = threshold + transitions.get_weight(context, i);
		sum += p[i];
	}

    real invsum = 1.0 / sum;
	for (int i=0; i<N; ++i) {
        p[i] *= invsum;
	}
}


/// Get the density of the transition probabilities \c p from \c context
real FactoredMarkovChain::pdf(Context context, Vector p)
{
	assert((context>=0)&&(context<tot_states));
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

