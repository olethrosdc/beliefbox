// -*- Mode: c++ -*-

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

template <class S, class A>
class Trajectory
{
public:
    std::vector<std::pair<S, A> > x;
    void Observe(S s, A a)
    {
        x.push_back(std::pair<S, A>(s, a));
    }
    
};

#endif
