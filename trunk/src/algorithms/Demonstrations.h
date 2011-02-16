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
    {
        NewEpisode();
    }
    void Observe(S s, A a)
    {
        current_trajectory->Observe(s, a);
    }
    void NewEpisode()
    {
        trajectories.push_back(Trajectory<S, A>());
        current_trajectory = &trajectories[trajectories.size()];
    }
};

#endif
