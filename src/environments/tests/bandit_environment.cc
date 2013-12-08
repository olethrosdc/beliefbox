/**
   Copyright (C) 2013, Christos Dimitrakakis
*/

#include <cstdio>
#include <cassert>
#include <string>
#include <vector>
#include <random>

// env_ function prototypes types 
#include <rlglue/Environment_common.h>	  

// helpful functions for allocating structs and cleaning them up 
#include <rlglue/utils/C/RLStruct_util.h> 

class BernoulliBandits
{
public:
	std::vector<double> means;
	std::default_random_engine generator;
	BernoulliBandits(int n) : means(n)
	{
		std::uniform_real_distribution<double> distribution(0.0,1.0);
		for (int i=0; i<n; ++i) {
			means[i] = distribution(generator);
		}
	}
	double GenerateReward(int a)
	{
		std::bernoulli_distribution(means[a]);
		return distribution(generator);
	}
};


// HELPER FUNCTIONS.  PROTOTYPES HERE, CODE AT THE BOTTOM OF FILE */

/* GLOBAL VARIABLES FOR RL-GLUE methods (global for convenience) */  
static world_description_t the_world;
static observation_t this_observation;
static reward_observation_terminal_t this_reward_observation;


/* Used if a message is sent to the environment to use fixed start states */

static string task_spec_string =  
	"VERSION RL-Glue-3.0 PROBLEMTYPE episodic \
DISCOUNTFACTOR 0.9999 OBSERVATIONS INTS (0 0) \
ACTIONS INTS (0 16)  REWARDS (0.0 1.0) \
EXTRA Bandits(C/C++) by Christos Dimitrakakis";


/** Initialise environment */

const char* env_init(){    
	the_world.numRows = 6;
	the_world.numCols = 18;

	/* Allocate the observation variable */
	allocateRLStruct(&this_observation, 1, 0, 0);
	/* That is equivalent to:
	   this_observation.numInts     =  1;
	   this_observation.intArray    = (int*)calloc(1,sizeof(int));
	   this_observation.numDoubles  = 0;
	   this_observation.doubleArray = 0;
	   this_observation.numChars    = 0;
	   this_observation.charArray   = 0;
	*/
	/* Setup the reward_observation variable */
	this_reward_observation.observation = &this_observation;
	this_reward_observation.reward=0;
	this_reward_observation.terminal=0;

	return task_spec_string.c_str(); 
}


/**
   Standard RL-Glue method. Sets an initial state and returns
   the corresponding observation.
*/
const observation_t *env_start()
{ 
	this_observation.intArray[0]=1;
  	return &this_observation;
}

const reward_observation_terminal_t *env_step(const action_t *this_action)
{
	/* Make sure the action is valid */
	assert(this_action->numInts==1);
	assert(this_action->intArray[0]>=0);
	assert(this_action->intArray[0]<4);

	updatePosition(&the_world,this_action->intArray[0]);
	this_reward_observation.observation->intArray[0] = calculate_flat_state(the_world);
	this_reward_observation.reward = calculate_reward(the_world);
	this_reward_observation.terminal = check_terminal(the_world.agentRow,the_world.agentCol);

	return &this_reward_observation;
}

void env_cleanup()
{
	clearRLStruct(&this_observation);
}

const char* env_message(const char* _inMessage) {
  
	string inMessage = _inMessage;

	if(inMessage == "test") {
        fixed_start_state=0;
        return "Test!";
    }
    

	/*	Message Description
 	 * 'print-state'
	 * Action: Print the map and the current agent location
	 */
	if(inMessage == "print-state"){
		print_state();
		return "Message understood.  Printed the state.";
	}

	return "SamplesMinesEnvironment(C++) does not respond to that message.";
}

