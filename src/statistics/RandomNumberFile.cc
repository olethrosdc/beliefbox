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

#include "RandomNumberFile.h"
#include <string>
#include <limits>
#include <cstdlib>
#include <cstdio>

RandomNumberFile::RandomNumberFile(const std::string & filename)
{
    FILE* fin = fopen (filename.c_str(), "r");
    if (!fin) {
        fprintf (stderr, "could not open %s for reading\n", filename.c_str());
        exit(-1);
    }
    
    
    fseek(fin, 0, SEEK_END);
    long flength = ftell(fin);
    rewind(fin);

    long length = sizeof(char)*flength/sizeof(ulong);
    pool.resize(length);

    fread((void*) &pool[0], sizeof(ulong), length, fin);

    fclose(fin);

    position = 0;
}


ulong RandomNumberFile::random()
{
    ulong x = pool[position];
    position = (position + 1) % pool.size();
    return x;
}

real RandomNumberFile::uniform() 
{
    const double LONG_INT_MAX = std::numeric_limits<ulong>::max();

    real x = 0;
    do {
        ulong y = random();
            //x = ((double) i) * (0.5 / (double) LONG_INT_MAX);
#if 0
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);
#endif
        x = ((real) y)/ LONG_INT_MAX;
//        printf("%f\n", x);
    } while (x >= 1.0);
    return x;
}
