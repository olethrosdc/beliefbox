/* -*- Mode: C++; -*- */
/* VER: $Id: ParticleFilter.h,v 1.3 2006/10/31 16:59:39 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include "real.h"
#include "Sampling.h"
#include "Distribution.h"
#include <vector>

class ParticleFilter
{
 public:
    int N;
    std::vector<real> y;
    std::vector<real> y2;
    std::vector<real> w;
    std::vector<real> w2;
    Distribution* transitions; ///< Transitions
    Distribution* observations; ///< Observations
    Distribution* prior; //<prior
    /// Constructor
	ParticleFilter(int N, Distribution* prior, Distribution* T, Distribution* O);
	void Init(int N, Distribution* prior, Distribution* T, Distribution* O);
    ParticleFilter();
    virtual void SetNumberOfEstimates(int n_estimates);
    virtual void Reset();
    virtual ~ParticleFilter();
    virtual real Sample();
    virtual void Observe(real x);
    virtual real GetMean();
    virtual real GetVar();
};


/// The problem is that ..
class BernoulliParticleFilter : public ParticleFilter 
{
 public:
    BernoulliParticleFilter(int N, Distribution* prior, Distribution* T, Distribution* O);
    virtual ~BernoulliParticleFilter();
    void Observe(real x);
    real GetMean();
    real GetVar();
};

class BernoulliGridParticleFilter : public ParticleFilter 
{
 public:
    /// Constructor
    BernoulliGridParticleFilter(int N, Distribution* prior, Distribution* T, Distribution* O);
    BernoulliGridParticleFilter();
    virtual ~BernoulliGridParticleFilter();
    void Observe(real x);
    real GetMean();
    real GetVar();
};



#endif
