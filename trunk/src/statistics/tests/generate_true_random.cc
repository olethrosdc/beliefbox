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
#include <cstdlib>

int main (int argc, char** argv)
{
    if (argc != 3) {
        fprintf(stderr, "Usage: generate_true_random numbers outputfile\n");
        return -1;
    }

    int numbers = atoi(argv[1]);
    char* outfname = argv[2];

    if (numbers <= 0) {
        fprintf (stderr, "numbers must be > 0\n");
        return -1;
    }

    FILE* fout = fopen(outfname, "w");

    if (!fout) {
        fprintf (stderr, "Could not open %s for writing\n", outfname);
        return -1;
    }


    for (int i=0; i<numbers; i++) {
        unsigned long x = true_random_bits();
        fwrite ((void*) &x, sizeof(unsigned long), 1, fout);
        printf("%f\n", 100.0 * (real) i / (real) numbers);
    }
    fclose (fout);
    return 0;
}

#endif
