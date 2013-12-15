/* -*- Mode: C++; -*- */
// copyright (c) 2013 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// -- Algorithms and models -- //
#include "SampleBasedRL.h"
#include "Random.h"
#include "MersenneTwister.h"
#include "DiscreteMDPCounts.h"

// -- Usual includes -- //
#include <cstring>
#include <getopt.h>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>

#include <rlglue/Agent_common.h> /* agent_ functions and RL-Glue types */
#include <rlglue/utils/C/RLStruct_util.h> /* helpful functions for structs */
#include <rlglue/utils/C/TaskSpec_Parser.h> /* task spec parser */



MersenneTwisterRNG rng;

class Agent
{
public:
  MDPModel* model;
  SampleBasedRL* sampling;
  action_t this_action;
  action_t last_action;
  observation_t *last_observation;
  int terminal_state;
  Agent(int n_states, int n_actions, real gamma)
    : last_observation(0),
      terminal_state(n_states)
  {
    double dirichlet_mass = 0.1;
    double epsilon = 0.0;
    int max_samples = 1;

    std::cout << "# Creating new MDP model"
              << std::endl;
    model =  new DiscreteMDPCounts(n_states + 1,
            n_actions,
            dirichlet_mass);

    std::cout << "# Creating new planner" << std::endl;
    sampling = new SampleBasedRL(n_states + 1,
            n_actions,
            gamma,
            epsilon,
            model,
            &rng,
            max_samples,
            true);

    std::cout << "# Allocating structures" << std::endl;
    allocateRLStruct(&this_action,1,0,0);
    allocateRLStruct(&last_action,1,0,0);
    last_observation=allocateRLStructPointer(1,0,0);
    std::cout << "# Initialisation complete" << std::endl;
  }
  ~Agent()
  {
    clearRLStruct(&this_action);
    clearRLStruct(&last_action);
    freeRLStructPointer(last_observation);

    delete sampling;
    delete model;
  }
  const action_t *Start(const observation_t *this_observation)
  {
    std::cout << "# Starting new episode" << std::endl;
    int state = this_observation->intArray[0];

    sampling->Reset();
    int theIntAction = sampling->Act(0.0, state);

    this_action.intArray[0] = theIntAction;
    replaceRLStruct(&this_action, &last_action);
    replaceRLStruct(this_observation, last_observation);
    return &this_action;
  }
  const action_t *Step(double reward, const observation_t *this_observation) {
    int newState = this_observation->intArray[0];
    //int lastState = last_observation->intArray[0];
    //int lastAction = last_action.intArray[0];
    int newAction = sampling->Act(reward, newState);

    this_action.intArray[0] = newAction;
    replaceRLStruct(&this_action, &last_action);
    replaceRLStruct(this_observation, last_observation);

    return &this_action;
  }
  void EpisodeEnd(double reward)
  {
    sampling->Act(reward, terminal_state);
  }
};

static Agent* agent;

/** Called at the beginning of an experiment */
void agent_init(const char* task_spec)
{
  /*Struct to hold the parsed task spec*/
  taskspec_t *ts= new taskspec_t; 
  int decode_result = decode_taskspec( ts, task_spec );
  if(decode_result!=0){
    std::cerr << "Could not decode task spec, code: " << decode_result
              << "for task spec: " << task_spec << std::endl; 
    exit(1);
  }
	
  // Lots of assertions to make sure that we can handle this problem.  
  assert(getNumIntObs(ts)==1);
  assert(getNumDoubleObs(ts)==0);
  assert(isIntObsMax_special(ts,0)==0);
  assert(isIntObsMin_special(ts,0)==0);
	
  int n_states =getIntObsMax(ts,0)+1;

  assert(getNumIntAct(ts)==1);
  assert(getNumDoubleAct(ts)==0);
  assert(isIntActMax_special(ts,0)==0);
  assert(isIntActMin_special(ts,0)==0);

  int n_actions = getIntActMax(ts,0)+1;

  double gamma = ts->discount_factor;

  free_taskspec_struct(ts); // Make the taskspec struct a "blank slate" 
  delete ts; // Free the structure itself 

  std::cout << "# Initialising new agent" << std::endl;

  agent = new Agent(n_states, n_actions, gamma);

}

/** Called at the beginning of an episode */
const action_t *agent_start(const observation_t *this_observation) {
  return agent->Start(this_observation);
}

/** Called at every step */
const action_t *agent_step(double reward, const observation_t *this_observation) {
  return agent->Step(reward, this_observation);
}

/** Called when the episode ends */
void agent_end(double reward) {
  return agent->EpisodeEnd(reward);
}

/** Called when the experiment ends */
void agent_cleanup()
{
  std::cout << "# Deleting agent" << std::endl;
  delete agent;
}

const char* agent_message(const char* _inMessage) 
{
  return "sample_agent_rl_glue does not understand your message.";

}
