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

/** Generic classifier class.
 */
template <typename InputSet, typename ClassSet, typename ClassDistribution>
class Classifier
{
public:
    /// Destructor
    virtual ~Classifier()
    {
    }
    /// Classify a specific input according to decision rule
    virtual ClassSet Classify(const InputSet& x) = 0;
    /// Return the complete class distribution for a specific input
    virtual ClassDistribution& Output(const InputSet& x) = 0;
    /// Observe a particular input and class label pair
    virtual real Observe(const InputSet& x, const ClassSet& y) = 0;
    /// Observe a particular input and class distribution
    virtual real Observe(const InputSet& x, const ClassDistribution& p) = 0;
};

#endif
