/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>

class ContextTreePredictor : public FactoredPredictor
{
protected:
    int n_actions; ///< the number of actions
    int n_obs; ///< the number of distinct observations
    ContextTree tree; ///< the context tree
    int current_obs; ///< the current observation
public:
    ContextTreePredictor(int n_actions_, int n_obs_, int depth)
        : n_actions(n_actions_),
          n_obs(n_obs_),
          tree(n_obs * n_actions, n_obs, depth),
          current_obs(0)
    {        
    }

    virtual ~ContextTreePredictor()
    {
    }
    /* Training and generation */
    /// Observe the (first?) observation.
    virtual real Observe (int prd) 
    {
        current_obs = prd;
        return 1.0 / (real) n_obs;
    }
    /// Observe current action and next observation
    virtual real Observe (int act, int prd) 
    {
        int x = act * n_obs + current_obs;
        current_obs = prd;
        return tree.Observe(x, prd);
    }

    virtual real ObservationProbability (int act, int x) 
    {
    }

    virtual void Reset()
    {
    }
};
