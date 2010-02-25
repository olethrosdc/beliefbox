// -*- Mode: c++ -*-
// copyright (c) 2010 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>

#ifndef RESOURCE_USE_H
#define RESOURCE_USE_H


#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

inline long getMaxResidentSetSize()
{
    static struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

inline long getSharedMemorySize()
{
    static struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_ixrss;
}

inline long getUnsharedDataSetSize()
{
    static struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_idrss;
}

inline long getUnsharedStackSize()
{
    static struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_isrss;
}

#endif
