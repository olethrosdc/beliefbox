/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef SPARSE_TRANSITIONS_H
#define SPARSE_TRANSITIONS_H

#include <map>
#include <vector>
#include "real.h"

/**
   \ingroup StatisticsGroup
 */
/*@{*/


typedef std::map<int, std::vector<real>, std::greater<int> > SourceMap;
typedef SourceMap::iterator SourceMapIterator;

/** A sparse transition model for discrete observations.
    
 */
class SparseTransitions 
{
protected:
	int n_sources; ///< number of source contexts
	int n_destinations; ///< number of next observations
	SourceMap sources; ///< a map of sources
public:
    /** Constructor

        @param n_sources number of contexts
        @param n_destinations number of predicted values
     */
	SparseTransitions(int n_sources, int n_destinations)
	{
		this->n_sources = n_sources;
		this->n_destinations = n_destinations;
	}

    /// Get the raw weight of a particular src/dst pair.
	real get_weight(int src, int dst) const
	{
		SourceMapIterator i = sources.find(src);
		if (i==sources.end()) {
			return 0.0;
		} else {
			return i->second[dst];
		}
	}
    /// Get the weights of all predictions from a src context.
    std::vector<real> get_weights(int src) const
    {
		SourceMapIterator i = sources.find(src);
		if (i==sources.end()) {
            std::vector<real> zero_vector(n_destinations);
			return zero_vector;
		} else {
			return i->second;
		}
	}
    /// Get the number of destinations
	int nof_destinations() const
	{
		return n_destinations;
	}

    /// Observe a particular transition
	real observe(int src, int dst)
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

/*@}*/
#endif
