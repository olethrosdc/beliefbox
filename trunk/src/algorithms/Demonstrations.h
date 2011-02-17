/// -*- Mode: c++ -*-

#ifndef DEMONSTRATIONS_H
#define DEMONSTRATIONS_H


#include "Trajectory.h"
template <class S, class A>
class Demonstrations
{
public:
    std::vector<Trajectory<S, A> > trajectories;
    Trajectory<S,A>* current_trajectory;
    Demonstrations() 
        : current_trajectory(NULL)
    {
        fprintf(stderr, "Creating Demonstrators\n");
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
};

#endif
