/* -*- Mode: C++; -*- */
/* VER: $Id: ParticleFilter.c,v 1.2 2006/11/06 18:28:10 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2004 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "ParticleFilter.h"
ParticleFilter::ParticleFilter(int N, Distribution* prior, Distribution* T, Distribution* O)
{
	this->transitions = T;
	this->observations = O;
	this->prior = prior;
	SetNumberOfEstimates (N);
}
    
void ParticleFilter::Init(int N, Distribution* prior, Distribution* T, Distribution* O)
{
	this->transitions = T;
	this->observations = O;
	this->prior = prior;
	SetNumberOfEstimates (N);
}

ParticleFilter::ParticleFilter()
{
	transitions = NULL;
	observations = NULL;
	prior = NULL;
	SetNumberOfEstimates(1);
}

void ParticleFilter::SetNumberOfEstimates(int n_estimates)
{
	N = n_estimates;
	y.resize(N);
	y2.resize(N);
	w.resize(N);
	w2.resize(N);
	if (prior) {
		Reset();
	}
}
    
void ParticleFilter::Reset()
{
    assert(prior);
	real a = 1.0f/(real) N;
	//printf ("Generating: ");
	for (int i=0; i<N; i++) {
		w[i] = a;
		y[i] = prior->generate();
		//printf ("(%f, %f)", w[i], y[i]);
	}
	//printf ("\n");
}

ParticleFilter::~ParticleFilter()
{
}

real ParticleFilter::Sample()
{
	int Yn = PropSample (w);
	return y[Yn] + transitions->generate(); 
}

void ParticleFilter::Observe(real x)
{
	// Generate a set of samples from our current belief
	for (int n=0; n<N; n++) {
		y2[n] = Sample();
	}

	// Evaluate the new posterior up to a normalising constant
	real sum = 0.0f;
	for (int i=0; i<N; i++) {
		real likelihood = observations->pdf(x - y2[i]);
		real prior = 0.0f;
		for (int j=0; j<N; j++) {
			prior += w[j]*transitions->pdf(y[j] - y2[i]);
		}
		w2[i] = likelihood*prior;
		sum += w2[i];
	}

	// Normalise to create a new filtering distribution
	real isum = 1.0f / sum;
	for (int i=0; i<N; i++) {
		w[i] = w2[i] * isum;
		y[i] = y2[i];
		//printf ("%f %f\n", w[i], y[i]);
	}
}

/// Get the current mean;
real ParticleFilter::GetMean()
{
	real mean = 0.0f;
	for (int i=0; i<N; i++) {
		mean += w[i] * y[i];
	}
	return mean;
}
    
/// Get the current mean;
real ParticleFilter::GetVar()
{
	real mean = GetMean();
	real var = 0.0f;
	for (int i=0; i<N; i++) {
		real delta = y[i] - mean;
		var += w[i] * delta * delta;
	}
	return var;
}


/// The problem is that ..
BernoulliParticleFilter::BernoulliParticleFilter(int N, Distribution* prior, Distribution* T, Distribution* O) : ParticleFilter(N, prior, T, O)
{
}

BernoulliParticleFilter::~BernoulliParticleFilter()
{
}

void BernoulliParticleFilter::Observe(real x)
{
	// Generate a set of samples from our current belief
	for (int n=0; n<N; n++) {
		int Yn = PropSample (w);
		y2[n] = y[Yn] + transitions->generate();
		w[n] = w[Yn];
		//printf ("y2:%f ", y2[n]);
	}
	
	// Evaluate the new posterior up to a normalising constant
	real sum = 0.0f;
	for (int i=0; i<N; i++) {
		real q = y2[i];
		real likelihood = q*x + (1-q)*(1-x);
		real prior = w[i];
		w2[i] = likelihood*prior;
		sum += w2[i];
	}
	if (sum==0.0f) {
		fprintf (stderr, "ERROR: 0 mass on prior!!\n");
		exit(-1);
	}

	// Normalise to create a new filtering distribution
	real isum = 1.0f / sum;
	for (int i=0; i<N; i++) {
		w[i] = w2[i] * isum;
		y[i] = y2[i];
		//printf ("%f %f\n", w[i], y[i]);
	}
}

/// Get the current mean;
real BernoulliParticleFilter::GetMean()
{
	real mean = 0.0f;
	for (int i=0; i<N; i++) {
		mean += w[i] * y[i];
	}
	return mean;
}
    
/// Get the current mean;
real BernoulliParticleFilter::GetVar()
{
	real mean = GetMean();
	real var = 0.0f;
	for (int i=0; i<N; i++) {
		real delta = y[i] - mean;
		var += w[i] * delta * delta;
	}
	return var;
}


BernoulliGridParticleFilter::BernoulliGridParticleFilter(int N, Distribution* prior, Distribution* T, Distribution* O) : ParticleFilter(N, prior, T, O)
{
}

BernoulliGridParticleFilter::BernoulliGridParticleFilter() : ParticleFilter()
{
}
BernoulliGridParticleFilter::~BernoulliGridParticleFilter()
{
}

void BernoulliGridParticleFilter::Observe(real x)
{
	// Evaluate the new posterior up to a normalising constant
	real sum = 0.0f;
	for (int i=0; i<N; i++) {
		real q = y[i];
		real likelihood = q*x + (1-q)*(1-x);
		real prior = w[i];
		w2[i] = likelihood*prior;
		sum += w2[i];
	}
	if (sum==0.0f) {
		fprintf (stderr, "ERROR: 0 mass on prior!!\n");
		exit(-1);
	}

	// Normalise to create a new filtering distribution
	real isum = 1.0f / sum;
	for (int i=0; i<N; i++) {
		w[i] = w2[i] * isum;
	}
}

/// Get the current mean;
real BernoulliGridParticleFilter::GetMean()
{
	real mean = 0.0f;
	for (int i=0; i<N; i++) {
		mean += w[i] * y[i];
	}
	return mean;
}
    
/// Get the current mean;
real BernoulliGridParticleFilter::GetVar()
{
	real mean = GetMean();
	real var = 0.0f;
	for (int i=0; i<N; i++) {
		real delta = y[i] - mean;
		var += w[i] * delta * delta;
	}
	return var;
}





