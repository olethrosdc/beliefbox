/* -*- Mode: C++; -*- */
/* VER: $Id: DiscretePolicy.h,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef DISCRETE_POLICY_H
#define DISCRETE_POLICY_H

#include "AbstractPolicy.h"


class DiscretePolicy : public AbstractPolicy<int, int>
{
public:
    DiscretePolicy() {}
    virtual ~DiscretePolicy() {}
    virtual int SelectAction() = 0;
    virtual void Observe (int previous_state, int action, real r, int next_state) = 0;
    virtual void Observe (real r, int next_state) = 0;
    virtual void Reset(int start_state) = 0;
    virtual real getActionProbability(int action) = 0;
    virtual real getActionProbability(int state, int action) = 0;
};

class FixedDiscretePolicy : public DiscretePolicy
{
public:
    std::vector<Vector> p;
    FixedDiscretePolicy (std::vector<Vector> p);
    virtual ~FixedDiscretePolicy();
    virtual int SelectAction();
    virtual void Observe (int previous_state, int action, real r, int next_state);
    virtual void Observe (real r, int next_state);
    virtual void Reset(int start_state = 0);
    virtual real getActionProbability(int action);
    virtual real getActionProbability(int state, int action);
    inline Vector getActionProbabilities(int state)
    {
        return p[state];
    }
    inline Vector* getActionProbabilitiesPtr(int state)
    {
        return &p[state];
    }
};

#endif
