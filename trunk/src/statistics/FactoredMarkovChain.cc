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
	  n_contexts((int) pow((double) n_states, (double) mem_size)),
      transitions(n_contexts, n_obs),
      act_history(mem_size),
      obs_history(mem_size),
      history(mem_size),
      threshold(0.5)
{
	
	printf("# Making FMC with %d actions, %d obs, %d states, %d history, %d contexts\n",
		   n_actions, n_obs, n_states, mem_size, n_contexts);
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
#if 0
    Context context = 0;
    Context n = 1;
    for (Context i=0; i<mem_size; i++, n*=n_states) {
        context += history[i]*n;
    }
#else
    Context context = history.get_id(n_states);
#endif
    return context;
}

/** Calculate the current context taking into account that there is an extra observation.

    Calculates the current contexts up to the time the
    last action was taken.
*/
FactoredMarkovChain::Context FactoredMarkovChain::getContext(int act)
{
    assert((obs_history.size() == 1 + act_history.size())
           && (act_history.size() == history.size()));
    // first set up the context component due to the most recent observation 
    // and postulated next action
	if (mem_size == 0) {
	  return 0;
	}
    Context context = CalculateState(act, obs_history.back());
    Context n = n_states;
    // continue with the remaining context
    for (Context i=0; i<mem_size - 1; i++, n*=n_states) {
        context += history[i]*n;
    }
    return context;
}


/** Obtain an observation.
    
    This should only be used at the first time step, before any actions are taken.

    @return Probability of having observed the next state.
*/
real FactoredMarkovChain::Observe(int prd)
{
	PushObservation(prd);
	return 0;
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
	//printf("%d %d %d %f # act obx ctx P\n", act, prd, current_context, Pr);
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
  assert((a>=0)&&(a<n_actions));
    
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
	assert((context>=0)&&(context<n_contexts));
    assert((prd>=0)&&(prd<n_states));
    return transitions.get_weight(context, prd);
}



/// Get the transition probability from \c context to \c prd
///
/// Takes into account the threshold.
real FactoredMarkovChain::getProbability(Context context, int prd)
{
    assert((context>=0)&&(context<n_contexts));
    assert((prd>=0)&&(prd<n_obs));
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
    assert((context>=0)&&(context<n_contexts));
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
        history[i] = 0;
        act_history[i] = 0;
        obs_history[i] = 0;
    }
    current_context = 0;
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



