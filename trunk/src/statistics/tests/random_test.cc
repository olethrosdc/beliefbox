/* -*- Mode: C++; -*- */
// copyright (c) 2008 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN
#include "Random.h"
#include "RandomNumberFile.h"
#include <cstdlib>

int main()
{
    fprintf (stderr,
             "sizes - real: %ld, float: %ld, double: %ld, int: %ld, ulong: %ld\n",
            sizeof(real),
            sizeof(float),
            sizeof(double),
            sizeof(int),
            sizeof(ulong));

    fprintf (stderr, "INT_MAX: %ld\n", INT_MAX);
    const ulong LONG_INT_MAX = (ulong) INT_MAX * (ulong) INT_MAX;
    fprintf (stderr, "LONG_INT_MAX: %ld\n", LONG_INT_MAX);
    fprintf (stderr, "ULONG_INT_MAX: %ld\n", 2*LONG_INT_MAX - 1);

    RandomNumberFile rng("/home/olethros/dev_random.bin");
    fprintf (stderr, "Rand int: %ld\n", rng.random());
    fprintf (stderr, "Rand real: %f\n", rng.uniform());
    //ulong max_long = 0;
    for (int i=0; i<rng.pool_size(); i++) {
        //ulong x = rng.random();
        //if (x > max_long) {
        //    max_long = x;
        // }
        printf ("%f\n", rng.uniform());
    }
    //fprintf (stderr, "%lu\n", max_long);
    return 0;
}


#endif




