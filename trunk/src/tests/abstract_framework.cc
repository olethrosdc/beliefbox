/* -*- Mode: C++; -*- */
/* VER: $Id: Sampling.h,v 1.3 2006/10/21 20:03:01 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
void DiscreteBelief::Update(Observation* observation)
{
	dynamic_cast<DiscreteObservation> observation
}
Action Policy::Act(Observation* observation)
{
	belief->Update(Observation* observation);
	return belief->Evaluate();
}

int main()
{
	action = policy->Act(observation);
	environment->PerformAction(action);
	
	return 0;
}
#endif
