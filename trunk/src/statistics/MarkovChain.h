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
#ifndef MARKOVCHAIN_H
#define MARKOVCHAIN_H

#include "Vector.h"
#include "Ring.h"
#include <vector>


/**
   \ingroup StatisticsGroup
 */
/*@{*/




#define EFFICIENT_MC_STATE_PUSH

/// A Markov Chain
class MarkovChain {
protected:
    typedef long MCState;
    int mem_pos; ///< memory position
	int n_states; ///< number of distinct states
	int mem_size; ///< history size for new transitions
	MCState curr_state; ///< current address in history
	MCState tot_states; ///< total number of representable states (mem_size*n_states)
#ifdef EFFICIENT_MC_STATE_PUSH
	Ring<int> memory; ///< hold history in here
#else
	std::vector<int> memory; ///< hold history in here
#endif
    MCState CalculateStateID ();
public:
    MarkovChain (int n_states, int mem_size);
    virtual ~MarkovChain ();

    /* probabilities */
    virtual real getTransition (MCState src, int dst) = 0;
    virtual real getProbability (MCState src, int dst) = 0;
    virtual void getProbabilities(MCState src, std::vector<real>& p) = 0;
    virtual void getNextStateProbabilities(std::vector<real>& p) = 0;
    virtual real pdf(MCState src, Vector q) = 0;
    virtual void setTransition (MCState src, int dst, real value) = 0;
    virtual void setThreshold (real threshold) = 0;

    /* Training and generation */
    int PushState (int state);    
    virtual real ObserveNextState (int state) = 0;
    virtual real NextStateProbability (int state) = 0;
    virtual void Reset () = 0;
    virtual int GenerateStatic () = 0;
    virtual int generate () = 0;
    MCState getCurrentState()
    {
        curr_state =  CalculateStateID ();
        return curr_state;
    }

    /* misc */
    MCState getTotalStates() {return tot_states;}
    /* Debug functions */
    virtual int ShowTransitions () = 0;
    void ShowMemory();
};




/*@}*/
#endif
