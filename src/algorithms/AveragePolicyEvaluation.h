// -*- Mode: c++ -*-
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef AVERAGE_POLICY_EVALUATION_H
#define AVERAGE_POLICY_EVALUATION_H

#include "PolicyEvaluation.h"

class AveragePolicyEvaluation : public PolicyEvaluation
{
public:
    AveragePolicyEvaluation(DiscretePolicy* policy_,
                            const DiscreteMDP* mdp_,
                            real baseline_ = 0.0);
    virtual ~AveragePolicyEvaluation();
    virtual void ComputeStateValues(real threshold, int max_iter=-1);
};

#endif

