/* -*- Mode: C++; -*- */
/* VER: $Id: MDPStateActionEstimator.c,v 1.1 2006/10/23 08:33:32 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "MDPStateActionEstimator.h"

MDPStateActionEstimator::MDPStateActionEstimator(MDPModel* model, real gamma, real init_val) 
    : StateActionEstimator(model->GetNStates(), model->GetNActions())
{
    this->model = model;
    //StateActionEstimator(model->GetNStates(), model->GetNActions());
}

MDPStateActionEstimator::~MDPStateActionEstimator()
{
}

void MDPStateActionEstimator::Reset()
{
}

void MDPStateActionEstimator::Observe(int s_p, int a_p, real r, int s, int a)
{
}

real MDPStateActionEstimator::TransProbability(int s_p, int a_p, int s)
{
    return model->getTransitionProbability(s_p, a_p, s);
}

int MDPStateActionEstimator::SampleTransition(int s, int a)
{
    return model->GenerateTransition(s, a);
}

real MDPStateActionEstimator::GetExpectedReward(int s, int a)
{
    return model->getExpectedReward(s, a);
}
real MDPStateActionEstimator::SampleReward(int s, int a)
{
    return model->getExpectedReward(s, a);
}

