// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BELIEF_DYNAMIC_PROGRAMMING_H
#define BELIEF_DYNAMIC_PROGRAMMING_H

/** Dynamic programming using beliefs.

This is dynamic programming using only beliefs, so no other state
variable is included.  

*/
template <typename Belief, typename Policy>
class BeliefDynamicProgramming
{
protected:
	Belief belief; ///< A belief is a certain type of distribution
	Policy policy;  
public:
	BeliefDynamicProgramming(Belief& initial_belief, Policy& policy)
		:
		belief(initial_belief)
		~virtual BeliefDynamicProgramming()
	{
	}
	void Expand(int T)
	{
		
	}
	
};

#endif

