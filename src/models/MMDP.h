// -*- Mode: c++ -*-
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: MDP.h,v 1.3 2006/11/06 23:42:32 olethros Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MMDP_H
#define MMDP_H

class Distribution;

#include "real.h"
#include "JointAction.h"


/** Abstract multi-agent MDP class */
class AbstractMMDP
{
	public:
    virtual ~AbstractMMDP() {}
};


/** Template for multi-agent MDP classes. */
template <typename StateType, typename ActionType>
class MMDP: public AbstractMMDP
{
 protected:
	StateType state;
 public:
	virtual ~MMDP()
	{
		// nothing to do
	}
	StateType getState() const
	{
		return state;
	}

};

#endif
