// -*- Mode: c++ -*-

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>

template <class S, class A>
class Trajectory
{
public:
    std::vector<S> x;
    
    Trajectory()
    {
    }

    void Observe(S s, A a)
    {
        //std::cerr << "Adding: " << s << ", " << a << std::endl;
        x.push_back(s);//std::pair<S, A>(s, a));
    }
    
};

#endif
