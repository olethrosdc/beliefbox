// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: MDP.h,v 1.3 2006/11/06 23:42:32 olethros Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef REWARD_DISTRIBUTION_H
#define REWARD_DISTRIBUTION_H

template <typename StateType, typename ActionType>
class RewardDistribution
{
 public:
    virtual ~RewardDistribution() {}
    virtual real generate(StateType state, ActionType action) = 0;
    virtual real expected(StateType state, ActionType action) = 0;
    virtual real pdf(StateType state, ActionType action, real reward) = 0;
};

#endif
