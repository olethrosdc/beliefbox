/* -*- Mode: C++; -*- */
// copyright (c) 2004-2011 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef WISHART_H
#define WISHART_H


#include "Distribution.h"

/// Wishart probability distribution
class Wishart : public VectorDistribution
{
public:
    Wishart();
    virtual ~Wishart();
    virtual void generate(T& x) const;
    virtual T generate() const;
    virtual real pdf(const T& x) const =0;
    virtual real log_pdf(const T& x) const
    {
        return log(pdf(x));
    }

};
