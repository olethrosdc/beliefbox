/*
 *  SinModel.cc
 *  
 *
 *  Created by Poseidon on 05/12/2013.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "SinModel.h"
#include "Random.h"
#include "RandomSourceRNG.h"
#include <cmath>

SinModel::Parameters SinModel::default_parameters = 
  {
    4,            // Upper bound on position
    -4,           // Lower bound on position
    0
  };

SinModel::SinModel(bool random_parameters)
  : Environment<Vector, int>(1, 1),
    parameters(default_parameters),
    noise(0.0, 0.1)
{
  if (random_parameters) {
    //RandomSourceRNG rng(false);
    parameters.U_POS = default_parameters.U_POS;
    parameters.L_POS = default_parameters.L_POS;
    parameters.MCNOISE = urandom() * default_parameters.MCNOISE;
  }
  state.Resize(n_states);
  state.Clear();
	
  state_upper_bound.Resize(n_states);
  state_lower_bound.Resize(n_states);
  state_upper_bound[0] = parameters.U_POS;
  state_lower_bound[0] = parameters.L_POS;
	
  action_upper_bound.Resize(n_actions);
  action_lower_bound.Resize(n_actions);
  for (uint i=0; i<n_actions; ++i) {
    action_lower_bound(i) = 0;
    action_upper_bound(i) = 0;
  }
    
  state_action_lower_bound.Resize(n_states + n_actions);
  state_action_upper_bound.Resize(n_states + n_actions);
  for (uint i=0; i<n_states; ++i) {
    state_action_lower_bound(i) = state_lower_bound(i);
    state_action_upper_bound(i) = state_upper_bound(i);
  }
  for (uint i=0; i<n_actions; ++i) {
    state_action_lower_bound(i + n_states) = action_lower_bound(i);
    state_action_upper_bound(i + n_states) = action_upper_bound(i);
  }
  endsim = false;
}

SinModel::~SinModel()
{
  // nothing to do
}


void SinModel::Reset()
{
  state[0] = urandom(parameters.L_POS, parameters.U_POS);
	
  endsim = false;
  reward = 0.0;
}
bool SinModel::Act(const int& action)
{
  // make sure we tell the guy we have terminated
  if (endsim) {
    reward = 0.0;
    return false;
  }
    
  // run
  Simulate(action);
	
  if (endsim) {
    return false;
  }
  return true;
}

void SinModel::Simulate(const int action)
{
  real input=0.0;
	
  switch (action){
  case 0: input = -1.0; break;
  default: Serror("Undefined action %d\n", action);
  }
    
  state[0] = sin(state[0]) + noise.generate();
    
  reward = 0;
  endsim = false;
    
  return;
	
}

