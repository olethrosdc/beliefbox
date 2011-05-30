// -*- Mode: c++ -*-

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>

template <class S, class A>
class Trajectory
{
public:
    std::vector<std::pair<S, A> > x;
    
    Trajectory()
    {
    }

    void Observe(S s, A a)
    {
        //std::cerr << "Adding: " << s << ", " << a << std::endl;
        x.push_back(std::pair<S, A>(s, a));
    }
    
	uint size() const
	{
		return x.size();
	}
	
	S state(uint t) const
	{
		assert (t < x.size());
		return x[t].first;
	}

	A action(uint t) const
	{
		assert (t < x.size());
		return x[t].second;
	}
};

#endif
