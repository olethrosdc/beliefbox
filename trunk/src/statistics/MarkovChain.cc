/* -*- Mode: c++ -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <stdexcept>
#include "debug.h"
#include "MarkovChain.h"
#include "Random.h"
#include "Ring.h"


/** \file MarkovChain.cc

	\brief Implementation of Markov chains.

	Interesting.
*/

/**
   \brief Creates a new markov chain.

   \arg \c n_states : The distinct observations possible. For example,
   if we wish to model transitions between letters, \c n_states might be
   equal to the number of letters in the alphabet, plus the number of
   punctuation symbols we use to incorporate.

   \arg \c mem_size : This is the Order of the Markov Chain. It
   specifies how large the state history should be. If 1, then only
   the previous observation is used to calculate transition
   probabilities. If n, then transition probabilities are considered
   between n-tuples of observation and the next observation.


   \return A pointer to the newly created chain, NULL if failed.
*/
MarkovChain::MarkovChain (int n_states, int mem_size)
{
    //Vector::BoundsCheckingStatus bounds_checking = Vector::NO_CHECK_BOUNDS;
    
    if (mem_size) {
        memory.resize(mem_size);
    }

	real max_bits = sizeof(MCState)*8;
	real used_bits = log((real) n_states) * ((real) (mem_size)) / log(2.0);
	
	if (max_bits < used_bits) {
		fprintf(stderr, "MC state is %f bits long, while %f bits are required for memory of length %d for %d states. Aborting\n", max_bits, used_bits, mem_size, n_states);
		exit(-1);
	} else {
		fprintf(stderr, "MC state is %f bits long, and %f bits are required for memory of length %d for %d states. Proceeding\n", max_bits, used_bits, mem_size, n_states);
	}
    mem_pos = 0;

	this->n_states = n_states;
	this->mem_size = mem_size;
	tot_states = 1;

	/* Clear up the memory before use! */
    memory.clear();

	for (int i=0; i<mem_size; i++) {
		tot_states *= n_states;
	}
    curr_state = 0;

}

//==========================================================
// MarkovChainCalculateStateID()
//----------------------------------------------------------
/**
   \brief Calculate a unique ID corresponding to the state history.

   Calculates a unique state ID from the MC's state history.
   The calculation is \f$\sum_i s_{t-i} N^i\f$, where \f$N\f$ is the number of
   states, \f$i \in [0,M-1]\f$, where \f$M\f$ is the history size (the order of
   the markov chain) and \f$s_{t}\f$ is the state at time \f$t\f$.

   \return The id of the current state history.

   SEE ALSO

   - MarkovChainPushState(), MarkovChainReset()
*/
MarkovChain::MCState MarkovChain::CalculateStateID () {
	MCState id = 0;
	MCState n = 1;
#ifdef EFFICIENT_MC_STATE_PUSH
    id = memory.get_id(n_states);
    //for (int i=0; i < mem_size; i++, n*=n_states) {
    //id += memory[i]*n; 
    //}
#else
	for (MCState i=1; i<=mem_size; i++, n*=n_states) {
		id += memory[i-1]*n; 
	}
#endif
	return id;
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
int  MarkovChain::PushState (int state) {
	int i;
	int popped_state;

    if (mem_size==0) { 
        return 0;
    }
#ifdef EFFICIENT_MC_STATE_PUSH
	popped_state = memory.front();
    memory.push_back(state);
#else
	popped_state = memory[mem_size - 1];
	for (i = mem_size - 1; i>0; i--) {
		memory[i] = memory[i-1];
	}
	memory[0] = state;
#endif    
    //logmsg("Pushing %d, popping %d\n", state, popped_state);
	return popped_state;
}

/**
   \brief Frees everything related to the markov chain.
*/
MarkovChain::~MarkovChain ()
{
}



void MarkovChain::ShowMemory()
{
    for (int i=0; i<mem_size; ++i) {
        printf ("%d ", memory[i]);
    }
    printf("\n");
}
