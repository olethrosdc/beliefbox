// -*- Mode: c++ -*-
// copyright (c) 2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef POLICY_BELIEF_H
#define POLICY_BELIEF_H



/** A belief about rewards in discrete spaces.
    
 */
class DiscretePolicyBelief
{
protected:
    int n_states;
    int n_actions;
public:
    // only set up
    DiscretePolicyBelief(int n_states_, int n_actions_)
        : n_states(n_states_), n_actions(n_actions_)
    {}
        
    virtual ~DiscretePolicyBelief()
    {}

    virtual real Update(int state, int action);
    virtual real CalculatePosterior(std::vecotr<int> states, std::vector<int> actions);
};




/** A belief about rewards in discrete spaces.
    
 */
class DirichletProductPolicyBelief
{
protected:
    std::vector<DirichletDistribution> P; ///< dirichlet distribution
public:
    // only set up
    DirichletProductPolicyBelief(int n_states_, int n_actions_)
        : DiscretePolicyBelief(n_states_, n_actions_)
    {
    }
        
    virtual ~DiscretePolicyBelief()
    {}

    virtual real Update(int state, int action);
    virtual real CalculatePosterior(std::vecotr<int> states, std::vector<int> actions);
};



#define
