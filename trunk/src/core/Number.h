// -*- Mode: c++ -*-
// copyright (c) 2005 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Id: Number.h,v 1.1 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef NUMBER_H
#define NUMBER_H

/** A general number class  
*/
class Number : public Object
{
public:
  virtual ~Number() {}
};

class Discrete : public Number
{
public:
  int value;
  virtual ~Discrete() {}
};

class Continuous : public Number
{
public:
  std::vector<float> value;
  virtual ~Continuous() {}
};

#endif
