// -*- Mode: c++ -*-
/* VER: $Id: EasyClock.h,v 1.1 2006/08/26 10:59:33 olethros Exp cdimitrakakis $*/
// copyright (c) 2006 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef EASY_CLOCK_H
#define EASY_CLOCK_H

#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

inline double GetCPU()
{
    static struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return (double) usage.ru_utime.tv_sec + ((double) usage.ru_utime.tv_usec)/1000000.0;
}

#endif
