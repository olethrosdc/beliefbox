// -*- Mode: c++ -*-
// copyright (c) 2021 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CRITIC_H
#define CRITIC_H

/// Abstract critic class
template <typename S, typename A>
class Critic
{
public:
	virtual ~Critic()
	{
	}
	/// observe \f$r_{t+1}, s_{t+1}, a_{t+1}\f$.
	virtual real Observe(real reward, const S& next_state, const A& next_action) = 0;
	/// Evaluate \f$s\f$ given the current history.
	virtual real getValue(const S& s) const = 0;
	/// Evaluate \f$s, a\f$ given the current history.
	virtual real getValue(const S& s, const A& a) const = 0;
	/// Reset the history.
	virtual void Reset() = 0;
};

#endif
