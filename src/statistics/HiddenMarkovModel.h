/* -*- Mode: c++;  -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef HIDDEN_MARKOV_MODEL_H
#define HIDDEN_MARKOV_MODEL_H


/**
   \ingroup StatisticsGroup
 */
/*@{*/


/// A generic hidden Markov model
template <typename S, typename X>
class HiddenMarkovModel
{
protected:
public:
    virtual ~HiddenMarkovModel ();
    virtual S getCurrentState();
    virtual X generate();
    virtual X generate_static();
};




/*@}*/
#endif
