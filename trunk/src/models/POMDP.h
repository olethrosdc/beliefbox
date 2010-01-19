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

#ifndef POMDP_H
#define POMDP_H

#include "real.h"
#include "SmartAssert.h"
#include "MDP.h"


template<typename ObservationType, typename StateType, typename ActionType>
class POMDP : public MDP<StateType, ActionType> {
public:
    virtual ~POMDP()
    {
    }
};

