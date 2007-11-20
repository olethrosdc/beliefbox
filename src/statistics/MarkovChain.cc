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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdexcept>
#include "debug.h"
#include "MarkovChain.h"
#include "Random.h"

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
    Vector::BoundsCheckingStatus bounds_checking = Vector::NO_CHECK_BOUNDS;
    
    if (mem_size) {
        memory.resize(mem_size);
    }

	this->n_states = n_states;
	this->mem_size = mem_size;
	tot_states = 1;

	/* Clear up the memory before use! */
	for (int i=0; i<mem_size; i++) {
		memory[i] = 0;
	}

	for (int i=0; i<mem_size; i++) {
		tot_states *= n_states;
	}
    curr_state = 0;

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
