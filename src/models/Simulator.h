// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef SIMULATOR_H
#define SIMULATOR_H

typedef <typename Model, typename Policy>
class Simulator
{
protected:
    Model& model;
    Policy& policy;
    real gamma;
    real last_total_reward;
    real last_discounted_reward;
public:
    Simulator(Model& model_, Policy& policy_, real gamma_) :
        model(model_),
        policy(policy_),
        gamma(gamma_)
    {
    }
    bool Simulate(int T)
    {
        policy.Reset(model.getState());
        for (int t=0; t<T; ++t) {
            A a = policy.SelectAction();
            bool success = model.Act(a);
            X x = model.getState();
            real r = model.getReward();
            policy.Observe(r, x);
            
        }
    }
};
#endif
