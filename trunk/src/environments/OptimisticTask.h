// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef OPTIMISTIC_TASK_H
#define OPTIMISTIC_TASK_H

#include "Environment.h"
#include "RandomNumberGenerator.h"
#include "DiscreteMDP.h"

class OptimisticTask : public DiscreteEnvironment
{
protected:
    real epsilon;
    real delta;
public:
    OptimisticTask(real epsilon_, real delta_);
    
    virtual ~OptimisticTask()
    {}
    
    virtual void Reset();
    virtual bool Act(int action);

    virtual const char* Name()
    {
        return "Optimistic Task";
    }

    virtual DiscreteMDP* getMDP() const;
};

#endif
