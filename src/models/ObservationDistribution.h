// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef OBSERVATION_DISTRIBUTION_H
#define OBSERVATION_DISTRIBUTION_H

template <typename StateType, typename ObservationType>
class ObservationDistribution
{
 public:
    virtual ~ObservationDistribution() {}
    virtual ObservationType generate(StateType state) = 0;
    virtual ObservationType expected(StateType state) = 0;
    virtual real pdf(StateType state, ObservationType observation) = 0;
};

#endif
