/// -*- Mode: c++ -*-

#ifndef DEMONSTRATIONS_H
#define DEMONSTRATIONS_H


#include "Trajectory.h"
#include "Environment.h"
#include "AbstractPolicy.h"
#include "Random.h"

template <class S, class A>
class Demonstrations
{
public:
    std::vector<Trajectory<S, A> > trajectories;
    Trajectory<S,A>* current_trajectory;
    Demonstrations() 
        : current_trajectory(NULL)
    {
        //fprintf(stderr, "Creating Demonstrators\n");
        NewEpisode();
    }
    void Observe(S s, A a)
    {
        //fprintf(stderr, "Size of trajectories: %d\n", trajectories.size())
        current_trajectory->Observe(s, a);
    }
    void NewEpisode()
    {
        //fprintf(stderr, "Adding Episode in Trajectories\n");
        trajectories.push_back(Trajectory<S, A>());
        current_trajectory = &trajectories[trajectories.size() - 1];
    }
	void Simulate(Environment<S, A>& environment, AbstractPolicy<S, A>& policy, real gamma, int horizon)
	{
		environment.Reset();
		policy.Reset();
		NewEpisode();
		bool running = true;
		real discount = 1.0;
		real total_reward = 0.0;
		real discounted_reward = 0.0;
		int t = 0;
		do {
			S state = environment->getState();
			real reward = environment->getReward();
			total_reward += reward;
			discounted_reward += discount * reward;
			discount *= gamma;
            A action = policy.Act(reward, state);
			Observe(state, action);
			running = environment->Act(action);
			if (horizon >= 0 && t >= horizon) {
				running = false;
			} else if (horizon < 0) {
				if (urandom() < 1.0 - gamma) {
					running = false;
				}
			}
		} while (running);
	}
    int size() const
    {
        if (trajectories[trajectories.size() - 1].size() > 0) {
            return trajectories.size();
        }
        return trajectories.size()  - 1;
    }
};

#endif
