/* -*- Mode: C++; -*- */
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "Vector.h"

int main(int argc, char** argv)
{
    
    {
        Vector x(4);
        for (uint i=0; i<x.Size(); ++i) {
        x(i) = i;
        }
        x.print(stdout);
        x *= 2;
        x.print(stdout);
        x = 3;
        x.print(stdout);
    }

    {
        real n = 3.0;
        Vector w(n);
        w.print(stdout);
    }

    {
        int n = 3;
        Vector w(n);
        w.print(stdout);

    }


    return 0;
}

#endif
