/* -*- Mode: c++;  -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "SparseMarkovChain.h"
#include "Random.h"

SparseMarkovChain::SparseMarkovChain(int n_states,
                                     int mem_size)
	: MarkovChain(n_states, mem_size),
      transitions((int) pow((double) n_states, (double) mem_size), n_states)
{

}

SparseMarkovChain::~SparseMarkovChain()
{
}

/**
   \brief Train a markov chain.

   - Incrementally trains the markov chain with a new
   observations. When wishing to add the statistics of a particular
   sequence to the transition probabilities, first call
   MarkovChainReset() and then call MarkovChainTrain() iteratively,
   presenting elements of the sequence in order.

   \arg \c chain: A pointer to the chain that you wish to train.
   \arg \c state: The state ID.

   \return Probability of having observed the next state.
*/
real SparseMarkovChain::ObserveNextState (int state)
{
	assert((state>=0)&&(state<n_states));
	curr_state = CalculateStateID();
	real Pr = getProbability(curr_state, state);
	transitions.observe(curr_state, state);
    PushState (state);
    return Pr;
}


/**
   \brief Likelihood of an observation for the markov chain.

   See what the probability of an observation is according to the
   current model.

   \arg \c state: The state ID.

   \return The probability of the transition to that state.
*/
real SparseMarkovChain::NextStateProbability (int state)
{
	assert((state>=0)&&(state<n_states));

	curr_state = CalculateStateID ();
	assert (curr_state>=0 && curr_state<tot_states);
	return getProbability(curr_state, state);
}


void SparseMarkovChain::getNextStateProbabilities(std::vector<real>& p)
{
    curr_state = CalculateStateID();
    return getProbabilities(curr_state, p);
}



/// Get the number of transitions from \c src to \c dst
real SparseMarkovChain::getTransition (MCState src, int dst)
{
	assert((src>=0)&&(src<n_states));
	assert((dst>=0)&&(dst<n_states));
	return transitions.get_weight(src, dst);
}



/// Get the transition probability from \c src to \c dst
///
/// Takes into account the threshold.
real SparseMarkovChain::getProbability(MCState src, int dst)
{
	assert((src>=0)&&(src<tot_states));
	assert((dst>=0)&&(dst<n_states));
	real sum = 0.0;
	int N = transitions.nof_destinations();
	
	for (int i=0; i<N; ++i) {
		sum += transitions.get_weight(src, i);
	}
	
	return (transitions.get_weight(src, dst) + threshold) / (sum + threshold * ((real) N));
}

/// Get the transition probabilities
///
/// Takes into account the threshold.
void SparseMarkovChain::getProbabilities(MCState src, std::vector<real>& p)
{
	assert((src>=0)&&(src<tot_states));
	assert((int) p.size()== n_states);
	real sum = 0.0;
	int N = transitions.nof_destinations();
	
	for (int i=0; i<N; ++i) {
        p[i] = threshold + transitions.get_weight(src, i);
		sum += p[i];
	}

    real invsum = 1.0 / sum;
	for (int i=0; i<N; ++i) {
        p[i] *= invsum;
	}
}


/// Get the density of the transition probabilities \c p from \c src
real SparseMarkovChain::pdf(MCState src, Vector p)
{
	assert((src>=0)&&(src<tot_states));
	Swarning("Not implemented\n");
	return 0.0;
}



/// Set the transition probability from \c src to \c dst
void SparseMarkovChain::setTransition (MCState src, int dst, real value)
{
	assert((src>=0)&&(src<tot_states));
	assert((dst>=0)&&(dst<n_states));
        //transitions[dst * tot_states + src] = value;
}

/**
   \brief Set the threshold specifying the Dirichlet prior.

   \arg \c threshold A threshold, added to all statistics. This way
   the probability of each transition will be \f$P(i|j)=\frac{N_{i,j}
   + \theta}{\sum_k N_{k,j} + \theta}\f$, there \f$\theta\f$ is the
   threshold.

   See also MarkovChainSoftmaxNormalize()
*/
void SparseMarkovChain::setThreshold (real threshold)
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
void SparseMarkovChain::Reset ()
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
int SparseMarkovChain::GenerateStatic ()
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

/**
   \brief Generate the next state
  
   Generates a value for the next state.

  
*/
int SparseMarkovChain::generate ()
{
	real tot = 0.0f;
    //int curr = CalculateStateID ();
	real sel = urandom();
	for (int j=0; j<n_states; j++) {
		real P = 0.0;//Pr[curr + j*tot_states];
		tot += P;
		if (sel<=tot && P>0.0f) {
			PushState(j);
			return j;
		}
	}
	int j = rand()%n_states;
	PushState(j);
	return j;
}




//==========================================================
// MarkovChainPushState()
//----------------------------------------------------------
/**
   Adds a new state to the history FIFO, pushing out the oldest
   member of the history. Note that the implementation could be
   faster, but the history is pretty small. If you change the
   implementation you also need to change
   MarkovChainCalculateStateID() at least.

   \arg \c chain: A pointer to the chain.
   \arg \c state: A new state

   \return It returns the ID of the popped state.
*/
int  SparseMarkovChain::PushState (int state) {
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

/**
   Can be useful for debugging.
   
   A way to test the Markov Chain's operation is to create a Markov
   Chain populated with randomly initialised transitions, use it to
   generate a sequence and then make a Markov Chain of the same order
   and train it on the generated data. The estimated transition
   probabilities should be close to those of the chain that generated
   the data.
*/

int SparseMarkovChain::ShowTransitions () {
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

