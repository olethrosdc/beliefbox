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
#include <iostream>
#include <cassert>
#include "SmartAssert.h"

// forward declaration
//template <typename C>
//class Ring;

/** A ring buffer.

    A ring buffer is useful for many things.
 */
template <typename C>
class Ring
{
public:
    class iterator
    {
    protected:
        const Ring<C>* ring;
        int index;
        int max_size;
    public:
        iterator()
        {
        }
        iterator(Ring<C>* ring_)
            : ring(ring_),
              index(ring->pos),
              max_size(ring->max_size)
        {
        }
        void operator--()
        {
            index++;
            while (index >= max_size) {
                index -= max_size;
            }
        }
        void operator++()
        {
            index--;
            while (index < 0) {
                index += max_size;
            }
        }
        int operator*()
        {
            assert(index >=0 && index < max_size);
            return ring->data[index];
        }
        bool operator==(const iterator& rhs) const
        {
            return (ring == rhs.ring && index == rhs.index);
        }
        bool operator!=(const iterator& rhs) const
        {
            return (ring != rhs.ring || index != rhs.index);
        }
           
    };
protected:
    C null_element;
    int max_size; ///< the maximum size of the array
    std::vector<C> data; ///< the data
    int pos; ///< the current position of the last write in the ring buffer. 
    ulong T;
    iterator end_iterator;
public:
    /// Construct a buffer of a given size
    Ring(int max_size_) : null_element(0),
                          max_size(max_size_),
                          data(max_size + 1),
                          pos(0),
                          T(0),
                          end_iterator(this)
    {
        clear();
    }
    
    void clear()
    {
        for (int i=0; i<data.size(); ++i) {
            data[i] = 0;
        }
        pos = 0;
        end_iterator.ring_pos = max_size;
    }
    /// Construct a buffer of a given size
    Ring() : null_element(0), max_size(0), pos(0), T(0)
    {
    }
    /// Destructor.
    ~Ring()
    {
    }
    iterator begin()
    {
        return (iterator (this));
    }
    iterator end()
    {
        return end_iterator
    }
    void resize(size_t new_size)
    {
        max_size = new_size;
        data.resize(max_size + 1);
        clear();
    }

    /** Get the unique ID of a state */
    long long get_id(int n_states)
    {
        if (max_size ==0) {
            return 0;
        }
        int n = 1;
        long long id = 0;
        for (int i = pos; i>=0; --i, n*=n_states) {
            id += data[i] * n;
        }
        
        for (int i = max_size - 1; i>pos; --i, n*=n_states) {
            id += data[i] * n; 
        }
        return id;
    }
    /** Push back a datum.
        
        Whenever new data is available, it is inserted in the current
        position. The position pointer wraps around.

        @param x The data to be pushed back.
     */
    void push_back(C x)
    {
        if (max_size) {
            pos++;
            if (pos >= (int) data.size()) {
                pos -= data.size();
            }
            data[pos] = x;
            int end_pos = pos + 1;
            if (end_pos > (int) data.size()) {
                end_pos = 
            }
        }
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
        int x = (pos - offset);
		/* Alternative implementations

		  1:  x - y*(x/y)
		  2: n = x % abs(y); if (n<0) n += abs(y);
          3: ((x % y) + y) & y
		 */
		x %= max_size;
		if (x < 0) 
		  x += max_size;

        //        std::cout << x << std::endl;
        SMART_ASSERT(x >= 0 && x < max_size)(x)(max_size);

        return data[x];
    }
    
    /// Return the oldest element in the ring buffer
    C& front()
    {
        if (max_size) {
            return data[(pos + 1)%max_size];
        } else {
            return null_element;
        }
    }

    /// Return the last element.
    C& back()
    {
        if (max_size) {
            return data[pos];
        } else { 
            return null_element;
        }
    }        
    ulong size()
    {
        return T;
    }
};

#endif
