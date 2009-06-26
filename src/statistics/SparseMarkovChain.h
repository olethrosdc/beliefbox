/* -*- Mode: c++;  -*- */
// copyright (c) 2007 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef SPARSE_MARKOVCHAIN_H
#define SPARSE_MARKOVCHAIN_H

#include "MarkovChain.h"
#include <map>
#include <vector>
#include <real.h>

/**
   \ingroup StatisticsGroup
 */
/*@{*/


typedef std::map<int, std::vector<real>, std::greater<int> > SourceMap;
typedef SourceMap::iterator SourceMapIterator;

class SparseTransitions 
{
protected:
	int n_sources;
	int n_destinations;
	SourceMap sources;
public:
	SparseTransitions(int n_sources, int n_destinations)
	{
		this->n_sources = n_sources;
		this->n_destinations = n_destinations;
	}
	real get_weight(int src, int dst)
	{
		SourceMapIterator i = sources.find(src);
		if (i==sources.end()) {
			return 0.0;
		} else {
			return i->second[dst];
		}
	}
    std::vector<real> get_weights(int src)
    {
		SourceMapIterator i = sources.find(src);
		if (i==sources.end()) {
            std::vector<real> zero_vector(n_destinations);
			return zero_vector;
		} else {
			return i->second;
		}
	}
	int nof_destinations()
	{
		return n_destinations;
	}

	float observe(int src, int dst)
	{
		SourceMapIterator i = sources.find(src);
		if (i==sources.end()) {
            std::vector<real> v(n_destinations);
            for (int j=0; j<n_destinations; ++j) {
                v[j] = 0.0;
            } 
            v[dst] = 1.0;
            std::pair<SourceMapIterator, bool> ret = sources.insert(std::make_pair(src, v));
            return 1.0;
        } else {
            i->second[dst]++;
            return i->second[dst];
        }
	}
};

/// A sparse implementation of a Markov chain
class SparseMarkovChain : public MarkovChain
{
protected:
    int n_transitions; ///< total number of transitions
	SparseTransitions transitions; ///< history-wide transition table

    float threshold;
public:
    SparseMarkovChain (int n_states, int mem_size);
    virtual ~SparseMarkovChain ();

    /* probabilities */
    virtual float getTransition (int src, int dst);
    virtual float getProbability (int src, int dst);
    virtual void getProbabilities(int src, std::vector<real>& p);
    virtual void getNextStateProbabilities(std::vector<real>& p);
    virtual float pdf(int src, Vector q);
    virtual void setTransition (int src, int dst, float value);
    virtual void setThreshold (float threshold);


    /* Training and generation */
    virtual int PushState (int state);
    virtual float ObserveNextState (int state);
    virtual float NextStateProbability (int state);
    virtual void Reset();
    virtual int GenerateStatic ();
    virtual int generate ();
    
    /* Debug functions */
    virtual int ShowTransitions ();
};
/*@}*/
#endif
