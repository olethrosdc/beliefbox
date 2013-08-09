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

#ifdef MAKE_MAIN
#include <vector>
#include "EasyClock.h"
#include "NormalDistribution.h"
#include "BetaDistribution.h"
#include "ProductDistribution.h"


int main (int argc, char** argv)
{
    NormalDistribution normal_distribution(0, 1.0);
    BetaDistribution beta_distribution(1.0, 1.0);


    std::vector<Distribution*> normal_beta(2);
    normal_beta[0] = &normal_distribution;
    normal_beta[1] = &beta_distribution;

    std::vector<Distribution*> normal_normal(2);
    normal_normal[0] = &normal_distribution;
    normal_normal[1] = &normal_distribution;

    std::vector<Distribution*> beta_beta(2);
    beta_beta[0] = &beta_distribution;
    beta_beta[1] = &beta_distribution;

    int T = 10;

    {
        printf(" # normal + beta\n");
        ProductDistribution<Distribution, real> product(normal_beta);
        for (int t=0; t<T; ++t) {
            Vector x(product.generate());
            std::vector<real> px = convert<real>(x);
            printf("%f | ", product.pdf(px));
            x.print(stdout);
            
        }
    }

    {
        printf(" # normal + normal\n");
        ProductDistribution<Distribution, real> product(normal_normal);
        for (int t=0; t<T; ++t) {
            Vector x(product.generate());
            std::vector<real> px = convert<real>(x);
            printf("%f | ", product.pdf(px));
            x.print(stdout);
        }
    }

    {
        printf(" # beta + beta\n");
        ProductDistribution<Distribution, real> product(beta_beta);
        for (int t=0; t<T; ++t) {
            Vector x(product.generate());
            std::vector<real> px = convert<real>(x);
            printf("%f | ", product.pdf(px));
            x.print(stdout);
        }
    }
    return 0;
}

#endif
