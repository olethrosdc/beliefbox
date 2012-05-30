// -*- Mode: c++ -*-
#ifndef ADVERSARY_H
#define ADVERSARY_H

#include "Matrix.h"
#include "DiscreteMDP.h"
#include "Dirichlet.h"

class Adversary
{
protected:
    Matrix counts;
    int n_states, n_actions;
    real gamma;
    Matrix reward;
    int last_state;
    int last_action;
    real last_reward;
public:
    Adversary(int n_states_, int n_actions_, real gamma_) :
        counts(n_states_, n_actions_), 
        n_states(n_states_),
        n_actions(n_actions_),
        gamma(gamma_),
        reward(n_states, n_actions),
        last_state(0),
        last_action(0),
        last_reward(0.0)
    {
    }
    virtual ~Adversary()
    {
    }
    void Observe(int state, int action)
    {
        counts(state, action)++;
        last_reward = reward(state, action);
        last_state = state;
        last_action = action;
    }
    virtual void setReward() = 0;
    real getReward(int state, int action) const
    {
        return reward(state, action);
    }
    real getReward() const
    {
        return last_reward;
    }
    const Matrix& getRewardMatrix() const
    {
        return reward;
    }
};

/** This adversary chooses a new, random reward, in every round */
class RandomAdversary : public Adversary
{
public:
    RandomAdversary(int n_states_, int n_actions_, real gamma_)
        : Adversary(n_states_, n_actions_, gamma_)
    {
    }
    virtual ~RandomAdversary()
    {
    }
    virtual void setReward()
    {
        DirichletDistribution Pr(n_states * n_actions, 1.0);
        Vector tmp_reward = Pr.generate();
        int i = 0;
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                reward(s, a) = tmp_reward(i++);
            }
        }	
    }
};


/** This adversary randomly chooses a fixed reward function in the beginning and never changes it. */
class FixedAdversary : public Adversary
{
public:
    FixedAdversary(int n_states_, int n_actions_, real gamma_)
        : Adversary(n_states_, n_actions_, gamma_)
    {
        DirichletDistribution Pr(n_states * n_actions, 1.0);
        Vector tmp_reward = Pr.generate();
        int i = 0;
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                reward(s, a) = tmp_reward(i++);
            }
        }	
    }
    virtual ~FixedAdversary()
    {
    }
    virtual void setReward()
    {
    }
};


class HeuristicAdversary : public Adversary
{
public:
    HeuristicAdversary(int n_states_, int n_actions_, real gamma_)
        : Adversary(n_states_, n_actions_, gamma_)
    {
    }
    virtual ~HeuristicAdversary()
    {
    }
    virtual void setReward()
    {
        //        real max_obs = Max(counts.RowMax());
        real min_obs = Min(counts.RowMin());
        //int i = 0;
        for (int s=0; s<n_states; ++s) {
            for (int a=0; a<n_actions; ++a) {
                if (counts(s,a) == min_obs) {
                    reward(s,a) = 1.0;
                } else {
                    reward(s,a) = 0.0;
                }
                //printf("%d %d %f %f #reward counts\n",
                //       s, a, reward(s,a), counts(s,a));
            }
        } 

    }
};






#endif
