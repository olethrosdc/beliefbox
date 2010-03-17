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

template <typename InputSet, typename ClassSet>
class AbstractClassifier
{
protected:
public:
    virtual ~AbstractClassifier() {}
    virtual ClassSet classify(InputSet x) = 0;
};

template <typename InputSet, typename ClassSet, typename Method>
class Classifier : public AbstractClassifier<InputSet, ClassSet>
{
protected:
    Method* method;
public:
    Classifier(Method* method_) : method(method_)
    {
    }
    virtual ~Classifier()
    {
    }
    virtual ClassSet classify(InputSet x)
    {
        return method->classify(x);
    }
};

#endif
