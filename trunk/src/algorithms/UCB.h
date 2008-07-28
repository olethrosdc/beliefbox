/* -*- Mode: C++; -*- */
/* VER: $Id: Policy.h,v 1.8 2006/10/23 08:33:24 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef UCB_POLICY_H
#define UCB_POLICY_H

#include <vector>
#include "real.h"
#include "DiscreteBanditPolicy.h"

class UCB1Policy : public DiscreteBanditPolicy
{
protected:
    std::vector<int> times;
    std::vector<real> Er;
    std::vector<real> U;
    int t;
    int n_actions;
public:
    UCB1Policy (int n_actions_);
    virtual ~UCB1Policy();
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual int SelectAction();
    void SetTimes(int a, int times_);
    void SetReward(int a, real reward);
};

class UCBgPolicy : public DiscreteBanditPolicy
{
protected:
    std::vector<int> times;
    std::vector<real> Er;
    std::vector<real> U;
    int t;
    int n_actions;
    real gamma;
public:
    UCBgPolicy (int n_actions_, real gamma_);
    virtual ~UCBgPolicy();
    virtual void Reset();
    virtual void Observe(int a, real r);
    virtual int SelectAction();
    void SetTimes(int a, int times_);
    void SetReward(int a, real reward);
};

class UCB2Policy : public DiscreteBanditPolicy
{
protected:
    std::vector<int> times;
    std::vector<real> Er;
    std::vector<real> U;
    int t;
    int n_actions;
public:
    UCB2Policy (int n_actions_);
    virtual ~UCB2Policy();
    virtual void Reset();
    virtual void Observe(int a, real r);

    virtual int SelectAction();
};



#endif
