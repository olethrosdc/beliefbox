// -*- Mode: c++ -*-
// copyright (c) 2013 by Nikolaos Tziortziotis <ntziorzi@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef SINMODEL_H
#define SINMODEL_H

#include "Environment.h"
#include "Vector.h"
#include "real.h"
#include "NormalDistribution.h"

/** Just produce a next from a sin + normal distribution
 */
class SinModel : public Environment<Vector, int>
{
protected:
  struct Parameters {
    real U_POS;         ///< Upper bound on position
    real L_POS;         ///< Lower bound on position
    real MCNOISE;       ///< input noise        
  };
  static Parameters default_parameters;
  Parameters parameters;
  Vector state_action_upper_bound;
  Vector state_action_lower_bound;
  Vector action_upper_bound;
  Vector action_lower_bound;
  void Simulate();
  NormalDistribution noise;
public:
  SinModel(bool random_parameters = false);
  virtual ~SinModel();
  virtual void Reset();
  virtual bool Act(const int& action);
  virtual void Simulate(const int action);
	
  const Vector& StateActionUpperBound() const
  {
    return state_action_upper_bound;
  }
  const Vector& StateActionLowerBound() const
  {
    return state_action_lower_bound;
  }
  const Vector& ActionUpperBound() const
  {
    return action_upper_bound;
  }
  const Vector& ActionLowerBound() const
  {
    return action_lower_bound;
  }
	
  virtual void setRandomness(real randomness)
  {
    parameters.MCNOISE = randomness;
  }
	
  virtual const char* Name() const
  {
    return "SinModel";
  }
	
  void Show()
  {
    printf("%f %f # params (SinModel)\n",
	   parameters.U_POS,
	   parameters.L_POS);
  }
  virtual real getTransitionProbability(const Vector& state, const int& action, const Vector& next_state) const
  {
    return 1.0;
  }
	
  virtual real getExpectedReward(const Vector& state, const int& action) const
  {
    return 0.0;
  }
};

class SinModelGenerator
{
public:
  SinModel Generate(bool random=true)
  {
    return SinModel(random);
  }
};







#endif
