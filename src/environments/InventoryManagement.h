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
#include "BinomialDistribution.h"
#include <string>
#include <vector>

// The inventory management task.
class InventoryManagement {
protected:
    int period;
    int max_items;
    real demand;
    std::vector<BinomialDistribution*> dist;
    DiscreteMDP* mdp;
public:
    InventoryManagement(int period_, int max_items_, real demand_);
    ~InventoryManagement();
    DiscreteMDP* getMDP() {return mdp;}
};

#endif
