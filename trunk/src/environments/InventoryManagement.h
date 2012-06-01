// -*- Mode: c++ -*-
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef INVENTORY_MANAGEMENT_H
#define INVENTORY_MANAGEMENT_H

#include "DiscreteMDP.h"
#include "Distribution.h"
#include "Environment.h"

#include <string>
#include <vector>


/// The inventory management task.
/// Described in Puterman, Chapter 3.
class InventoryManagement : public DiscreteEnvironment {
protected:
    int period; ///< This is the period (number of steps) between orders
    int max_items; ///< This is the maximum number of items that can be stocked
    real demand; ///< This is the probability that an item will be bought at each step
    real margin; ///< this is the ratio (greater than 1) of the retail value of the item versus the order cost
    //std::vector<SingularDistribution*> dist; ///< container for the reward distributions
    DiscreteMDP* local_mdp; ///< points to the MDP

public:
    InventoryManagement(int period_, int max_items_, real demand_, real margin_);
    ~InventoryManagement();
    DiscreteMDP* getMDP() const;
   
    virtual void Reset() 
    {
        state = 1;
        local_mdp->setState(state);
    }

    virtual bool Act(const int action) 
    {
        reward = local_mdp->Act(action);
        state = local_mdp->getState();
        return true;
    }
    virtual const char* Name()
    {
        return "Inventory management";
    }

};

#endif
