/* -*- Mode: C++; -*- */
/* VER: $Id: Distribution.h,v 1.3 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $*/
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/** An abstract type of prior distribution.

    Model a family of distributions on X, with the parameter of each
    distribution taking values in Q.

 */
template <typename X, typename Q>
class PriorDistribution
{
public:
    /// Update the parameter distribution, return the marginal probability
    virtual real Observe(const X& x) = 0;
    /// Return the marginal probability
    virtual real marginal(const X& x) const = 0;
    /// Return the probability of the particular family member
    virtual real pdf(const Q& q) const = 0;
    /// generate a member (or parameters) of the family
    virtual Q generate() = 0;
    /// generate observations
    virtual X generate_observations() = 0;
};



