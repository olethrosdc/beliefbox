/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CLASSIFIER_POLICY_H
#define CLASSIFIER_POLICY_H

#include "AbstractPolicy.h"

class ClassifierPolicy  : public AbstractPolicy<Vector, int>
{
protected:
    Classifier<Vector,int,Vector>* classifier;
public:
    ClassifierPolicy(Classifier<Vector,int,Vector>* classifier_) :
        classifier(classifier_)
    {
    }
    virtual ~ClassifierPolicy()
    {
    }
    virtual int SelectAction()
    {
        return classifier->Classify(state);
    }
    virtual void Observe (const Vector& previous_state, const int& action, real r, const Vector& next_state)
    {
        classifier->Observe(previous_state, action);
        this->state = next_state;
    }
    virtual void Observe (real r,  const Vector& next_state) 
    {
        state = next_state;
        return;
    }
    virtual void Reset()
    {
    }
    virtual void SetState(const Vector& state)
    { 
        this->state = state;
    }
};
#endif
