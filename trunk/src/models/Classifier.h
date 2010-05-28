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

template <typename InputSet, typename ClassSet, typename ClassDistriution>
class AbstractClassifier
{
protected:
public:
    virtual ~AbstractClassifier() {}
    virtual ClassSet Classify(InputSet& x) = 0;
};

template <typename InputSet, typename ClassSet, typename ClassDistriution, typename Method>
class Classifier : public AbstractClassifier<InputSet, ClassSet, ClassDistriution>
{
protected:
    Method* method;
public:
    Classifier(Method* method_) : method(method_)
    {
    }
    virtual ~Classifier
    {
    }
    virtual ClassSet Classify(InputSet& x)
    {
        return method->Classify(x);
    }
    virtual ClassDistribution& Output(InputSet& x)
    {
        return method->Output(x);
    }
    virtual void Observe(InputSet& x)
    {
        method->Observe(x);
    }
};

#endif
