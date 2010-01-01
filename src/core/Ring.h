/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RING_H
#define RING_H

#include <vector>

/** A ring buffer.

    A ring buffer is useful for many things.
 */
template <typename C>
class Ring
{
public:
    int max_size;
    std::vector<C> data;
    int pos;
    ulong T;
    /// Construct a buffer of a given size
    Ring(int max_size_) : max_size(max_size_), data(max_size), pos(0), T(0)
    {
        assert(max_size > 0);
        for (int i=0; i<max_size; ++i) {
            data[i] = 0;
        }
    }
    /// Destructor.
    ~Ring()
    {
    }
    /** Push back a datum.
        
        Whenever new data is available, it is inserted in the current
        position. The position pointer wraps around.

        @param x The data to be pushed back.
     */
    void push_back(C x)
    {
        pos = (pos + 1) % max_size;
        data[pos] = x;
        T++;
    }
    
    /** Access data at offset from the current position.
        
        If we are currently at time T, then we access the data
        inserted at time T - offset.
        
        @param offset This must be non-negative and smaller than \c max_size.

        @return A reference to the data.
     */
    C& operator[](int offset)
    {
        assert(offset >= 0 && offset < max_size);
        return data[(offset - pos) % max_size];
    }
    
    /// Return the last element.
    C& back()
    {
        return data[0];
    }

    int size()
    {
        return T;
    }
};

#endif
