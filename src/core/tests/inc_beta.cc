/* -*- Mode: c++ -*- */
/* VER: $Id: MathFunctions.h,v 1.2 2006/11/06 15:48:53 cdimitrakakis Exp cdimitrakakis $ */
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef MAKE_MAIN

#include "SpecialFunctions.h"
#include "EasyClock.h"
#include <cstdlib>
#include <cstdio>
#include <exception>
#include <stdexcept>

int main(void)
{
    
    float x = 0.9;
    float a = 7.0;
    float b = 0.5;
    printf ("betainc (%f, %f, %f) = %f\n", x, a, b, betai(x, a, b));
    return 0;
}

#endif
