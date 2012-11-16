// -*- Mode: c++ -*-

#ifndef MDP_WRAPPER_H
#define MDP_WRAPPER_H

#include "DiscreteMDP.h"
#include "Environment.h"

class MDPWrapper : public DiscreteEnvironment
{
protected:
  DiscreteMDP* mdp;
public:
  MDPWrapper(DiscreteMDP* mdp_) :
    Environment<int,int>(mdp_->getNStates(), mdp_->getNActions()),
    mdp(mdp_)
  {
    
  }
  virtual ~MDPWrapper()
  {
  }

  virtual void Reset()
  {
    state = rand()%n_states;
    reward = 0.0;
    mdp->setState(state);
  }
  virtual void Reset(int new_state)
  {
    state = new_state;
    reward = 0.0;
    mdp->setState(new_state);
  }

  virtual bool Act(const int action)
  {
    reward = mdp->Act(action);
    state = mdp->getState();
    return true;
  }

  virtual const char* Name()
  {
    return "Wrapper for some MDP";
  }
  
};

#endif
