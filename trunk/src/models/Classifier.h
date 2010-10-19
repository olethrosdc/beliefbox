/* -*- Mode: c++;  -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef CLASSIFIER_H
#define CLASSIFIER_H

template <typename InputSet, typename ClassSet, typename ClassDistribution>
class Classifier
{
public:
    virtual ~Classifier()
    {
    }
    virtual ClassSet Classify(const InputSet& x) = 0;
    virtual ClassDistribution& Output(const InputSet& x) = 0;
    virtual real Observe(const InputSet& x, const ClassSet& y) = 0;
};

#endif
