/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef FACTORED_PREDICTOR_H
#define FACTORED_PREDICTOR_H

/// Abstract class for prediction with actios
class FactoredPredictor
{
public:
  virtual ~FactoredPredictor()
  {}
    
  /* Training and generation */
  virtual real Observe (int prd) = 0;
  virtual real Observe (int act, int prd) = 0;
  virtual real ObservationProbability (int act, int x) = 0;
  //virtual real ObservationProbability (int x) = 0;
  virtual void Reset() = 0;
    
}; 



#endif
