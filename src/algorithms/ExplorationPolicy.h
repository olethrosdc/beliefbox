// -*- Mode: c++ -*-
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
// $Revision$
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef EXPLORATION_POLICY_H
#define EXPLORATION_POLICY_H

#include <cassert>
#include "Matrix.h"
#include "Random.h"
#include "Distribution.h"

/// Value-function-based exploration policy
///
/// Examples: epsilon-greedy, softmax
class VFExplorationPolicy
{
public:
    virtual ~VFExplorationPolicy()
    {}
    virtual void Observe(real reward, int state) = 0;
    virtual int SelectAction() = 0;
    virtual void setValueMatrix(const Matrix* Q) = 0;
    virtual DiscretePolicy* getFixedPolicy() = 0;
    virtual void setGeometricSchedule(real alpha, real beta)
    {}
};

class EpsilonGreedy : public VFExplorationPolicy
{
protected:
    int n_actions;
    real epsilon;
    int state;
    const Matrix* Q;
    bool use_geometric_schedule;
    real alpha, beta;
public:
    EpsilonGreedy(int n_actions_, real epsilon_) :
        n_actions(n_actions_), epsilon(epsilon_), Q(NULL),
        use_geometric_schedule(false)
    {
        assert(n_actions > 0);
        assert(epsilon >= 0 && epsilon <= 1);
    }

    virtual ~EpsilonGreedy()
    {
    }

    real getEpsilon()
    {
        return epsilon;
    }
    real setEpsilon(real epsilon_)
    {
        epsilon = epsilon_;
        assert(epsilon >= 0 && epsilon <= 1);
        return epsilon;
    }
    virtual void setValueMatrix(const Matrix* Q_)
    {
        Q = Q_;
    }

    virtual void setGeometricSchedule(real alpha_, real beta_)
        
    {
        alpha = alpha_;
        beta = beta_;
        use_geometric_schedule = true;
    }
    virtual int SelectAction() 
    {
        real threshold = epsilon;
        if (use_geometric_schedule) {
            threshold = epsilon / (1 + sqrt(beta));
            beta += alpha;
        }
        if (urandom() < threshold) {
            int action =  (int) floor(urandom(0.0, n_actions));
            //printf("returning random action %d\n", action);
            return action;
        }

        int argmax = 0;
        real max = (*Q)(state, argmax);
        for (int a=1; a<n_actions; ++a) {
            real Qsa = (*Q)(state, a);
            if (Qsa > max) {
                max = Qsa;
                argmax = a;
            }
        }
        return argmax;
    }
    virtual void Observe(real reward, int state)
    {
        this->state = state;
    }

    virtual DiscretePolicy* getFixedPolicy() 
    {
        Serror ("Not implemented\n");
        exit(-1);
        return NULL;
    }

};

class SoftmaxPolicy : public VFExplorationPolicy
{
protected:
    int n_actions;
    real beta;
    int state;
    const Matrix* Q;
public:
    SoftmaxPolicy(int n_actions_, real beta_) :
        n_actions(n_actions_), beta(beta_), Q(NULL)
    {
        assert(n_actions > 0);
        assert(beta >= 0);
    }

    virtual ~SoftmaxPolicy()
    {
    }

    real getBeta()
    {
        return beta;
    }
    real setBeta(real beta_)
    {
        beta = beta_;
        assert(beta >= 0);
        return beta;
    }
    virtual void setValueMatrix(const Matrix* Q_)
    {
        Q = Q_;
    }

    virtual int SelectAction() 
    {
		Vector eQ(n_actions);
		real s = 0;
		for (int a=0; a<n_actions; ++a) {
			eQ(a) = beta * exp((*Q)(state, a));
			s += eQ(a);
		}
		eQ /= s;
		return DiscreteDistribution::generate(eQ);
	}

    virtual void Observe(real reward, int state)
    {
        this->state = state;
    }

    virtual DiscretePolicy* getFixedPolicy() 
    {
        Serror ("Not implemented\n");
        exit(-1);
        return NULL;
    }

};

class MaxSoftmaxPolicy : public VFExplorationPolicy
{
protected:
    int n_actions;
    real beta;
	real epsilon;
    int state;
    const Matrix* Q;
public:
    MaxSoftmaxPolicy(int n_actions_, real beta_, real epsilon_) :
        n_actions(n_actions_), beta(beta_), epsilon(epsilon_), Q(NULL)
    {
        assert(n_actions > 0);
        assert(beta >= 0);
    }

    virtual ~MaxSoftmaxPolicy()
    {
    }

    real getBeta()
    {
        return beta;
    }
    real setBeta(real beta_)
    {
        beta = beta_;
        assert(beta >= 0);
        return beta;
    }
    real getEpsilon()
    {
        return epsilon;
    }
    real setEpsilon(real epsilon_)
    {
		epsilon = epsilon_;
        assert(epsilon >= 0 && epsilon <=1);
        return epsilon;
    }
    virtual void setValueMatrix(const Matrix* Q_)
    {
        Q = Q_;
    }

    virtual int SelectAction() 
    {
		// If we must select a random action, use Boltzmann instead of uniform.
		if (urandom() <= epsilon) {
			Vector eQ(n_actions);
			real s = 0;
			for (int a=0; a<n_actions; ++a) {
				eQ(a) = beta * exp((*Q)(state, a));
				s += eQ(a);
			}
			eQ /= s;
			return DiscreteDistribution::generate(eQ);
		}
		// If we must select a non-random action, select the maximum.
		int argmax = 0;
		real max = (*Q)(state, argmax);
		for (int a=1; a<n_actions; ++a) {
			real Qsa = (*Q)(state, a);
			if (Qsa > max) {
				max = Qsa;
				argmax = a;
			}
		}
		return argmax;
	}

    virtual void Observe(real reward, int state)
    {
        this->state = state;
    }

    virtual DiscretePolicy* getFixedPolicy() 
    {
        Serror ("Not implemented\n");
        exit(-1);
        return NULL;
    }

};
#endif
