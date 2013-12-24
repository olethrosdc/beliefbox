// -*- Mode: c++ -*-
// Various utility functions
#ifndef HASH_COMBINE_H
#define HASH_COMBINE_H

#include <unordered_map>

#define HASH_COMBINE 0x9e3779b9

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
	std::hash<T> hasher;
	seed ^= hasher(v) + HASH_COMBINE + (seed << 6) + (seed >> 2);
}


#endif
