// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ONLINE_ALGORITHM_H
#define ONLINE_ALGORITHM_H

template <typename A, typename S>
class OnlineAlgorithm
{
public:
    OnlineAlgorithm()
    {
    }
    virtual ~OnlineAlgorithm()
    {
    }
    /// call this at the end of an episode.
    virtual void Reset()
    {
    }
    /// Partial SARSA observation (can be used with eligibility traces)
    virtual real Observe (real reward, S next_state, A next_action) = 0;
    /// Get an action using the current exploration policy.
    /// it calls Observe as a side-effect.
    virtual A Act(real reward, S next_state) = 0;
    virtual real getValue (S state, A action) = 0;
};

#endif
