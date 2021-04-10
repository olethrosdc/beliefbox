/* -*- Mode: C++; -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <vector>
#include <cassert>
#include "Vector.h"
#include "Grid.h"
#include "DirichletTransitions.h"
#include "Matrix.h"

/** A simple radial basis function */
class RBF
{
public:
    Vector center;		///< The centroid
    Vector beta;		///< the variance
    /// Constructor
    RBF(const Vector& c, real b) : center(c)
    {
        assert(b > 0);
        beta.Resize(c.Size());
		for (int i=0; i<beta.Size(); i++) {
			beta(i) = b;
		}
    }
	
    RBF(const Vector& c, const Vector& b) : center(c), beta(b)
    {
        assert(b > 0);
    }
	
    /// Get the density at point x
    real Evaluate(const Vector& x)
    {
        Vector d = pow((x - center)/beta, 2.0);
        real r = d.Sum();
        return exp(-0.5*r);
    }
    /// Evaluate the log density
    real logEvaluate(const Vector& x)
    {
        Vector d = pow((x - center)/beta, 2.0);

        return d.Sum();
        //    return (-beta) * EuclideanNorm(&x, &center);		
    }
};

/// Abstract class for a basis
template <typename S, typename A>
class BasisSet
{
protected:
	int n_bases;
	mutable Vector log_features;
    mutable Vector features;
    mutable bool valid_features;
    mutable bool valid_log_features;
public:
	BasisSet()
		: n_bases(0),
		  valid_features(false),
		  valid_log_features(false)
	{
	}
	virtual ~BasisSet()
	{
	}
	int size()
    {
        return n_bases;
    }
	/// Reset statistics for history-dependent features
	virtual void Reset()  = 0;
	/// Obtain a new observation \f$a_t, r_t, s_{t+1}\f$.
	///
	/// To be used especially with history-dependent features
	virtual void Observe(const A& action, real reward, const S& next_state) = 0;
	/// Obtain state observation \f$s_1\f$.
	virtual void Observe(const S& state) = 0;
	real log_F(int j) const
    {
        assert(j >= 0 && j < n_bases);
        assert(valid_log_features);
        return log_features[j];
    }
    Vector log_F() const
    {
		assert(valid_log_features);
        return log_features;
    }
    real F(int j) const
    {
        assert(j >= 0 && j < n_bases);
        assert(valid_features);
        return features[j];
    }
    Vector F() const
    {
		assert(valid_features);
        return features;
    }
};

/// Simple RBF basis
///
/// Provides a basis only over states, ignoring actions.
class RBFBasisSet : public BasisSet<Vector, int>
{
protected:
    std::vector<RBF*> centers;
public:
    RBFBasisSet()
    {
		Reset();
	}
    RBFBasisSet(const EvenGrid& grid, real scale = 1);
    virtual ~RBFBasisSet();
	virtual void Reset()
	{
		valid_log_features = false;
		valid_features = false;
	}
	virtual void Observe(const int& action, real reward, const Vector& next_state);
	virtual void Observe(const Vector& state);
	
    void AddCenter(const Vector& v, const Vector& b);
    void AddCenter(const Vector& v, real b);
    void Evaluate(const Vector& x) const;
    void logEvaluate(const Vector& x) const;


};


/// Simple basis for counts
class CountsBasis : public BasisSet<int, int>
{
protected:
	int n_states;
	int n_actions;
	int state;
	DirichletTransitions model;
	Matrix average_reward;

public:
    CountsBasis(int n_states_, int n_actions_): 
		n_states(n_states_),
		n_actions(n_actions_),
		model(n_states, n_actions),
		average_reward(n_states, n_actions)
    {
		n_bases = n_states*n_actions*n_states * 2 + n_states*n_actions + n_states + 1;
		features = Vector(n_bases, Vector::CHECK_BOUNDS);
		Reset();
	}
	virtual ~CountsBasis();
	virtual void Reset();
	virtual void Observe(const int& action, real reward, const int& next_state);
};

#endif
